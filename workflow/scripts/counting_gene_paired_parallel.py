#!/usr/bin/env python3

##############################################################
# counting_gene_paired_parallel.py
# Description: This script attaches ligation and RT barcodes and UMI to FASTQ files
#   Input:
#      --input_folder:   Folder with BAM files
#      --output_folder:  Output folder for count files
#      --i7_prefix_file: Text file with i7 demultiplexed prefix
#      --gtf:            GTF file with exon information only
#      --rt_barcode_tsv: RT barcode tsv (Sample	shortdt	randomN)
#      --threads:        Number of threads
#      --gzip:           Compress the output files
#   Output:
#      <output_folder>/<i7_prefix>.count: Count file with gene expression per cell ┐ 
#      <output_folder>/<i7_prefix>_cell_id.csv: Cell id file with unmatched rate   │ for later sparse matrix creation
#      <output_folder>/gene_ids.csv: Gene exon id table                            ┘
##############################################################

# Counting Strategy:
# Get the intersection and union of genes from gene/exon arrays
#    1. If intersection is 1, count it
#    2. If intersection is more than 1, find the gene with the nearest transcript end (within 100bp)
#    3. If intersection is 0, try union
#    4. If union is 1, count it
#    5. If union is more than 1, find the gene with the nearest transcript end (within 100bp)
# 
# 1. Try to find by exons
# 2. Try to find by genes
# 3. Try to find by flanking regions (genes) iv - 1000 bp on + strand, iv + 1000 bp on - strand
# 4. Try to find on the opposite strand by exons
# 5. Try to find on the opposite strand by genes
# 6. If still not found, count as no feature

import collections
from datetime import datetime
import HTSeq
import os
import pandas as pd
import numpy as np
import gzip
IS_GZIP = False

def count_pair(algn_pair, genes, exons, transcript_ends, is_randomn, counter):
    """
    Count gene reads in a pair of alignments
    """
    if not algn_pair[0].aligned or not algn_pair[1].aligned:
        counter["_unmapped"] += 1 # Just in case, since we only have primary alignments
        return
        # Because Counter is a reference here, we don't need to return it.
    if algn_pair[0].iv.chrom not in genes.chrom_vectors or algn_pair[1].iv.chrom not in genes.chrom_vectors:
        counter["_unmapped"] += 1 # Just in case, since we use the same GTF file for alignment and counting
        return

    # Because of the Illumina sequencing and library preparation.
    #   1st read is always anti-sense to the gene(mRNA), while 2nd read is sense.
    reverse_strand(algn_pair[0])
    
    # Exon Match
    gene_ids_union, gene_ids_intersect = get_union_intersect_genes(algn_pair, exons)
    if len(gene_ids_intersect) == 1:
        counter[next(iter(gene_ids_intersect))] += 1
        return
    elif len(gene_ids_intersect) > 1:
        gene_id = find_nearest_transcript_end(algn_pair[0].iv.end_d, gene_ids_intersect, transcript_ends, is_randomn)
        counter[gene_id] += 1
        return
    
    # Intron Match
    suffix = "_intron" 
    if len(gene_ids_union) == 1:
        counter[next(iter(gene_ids_union))] += 1
        return
    elif len(gene_ids_union) > 1:
        gene_id = find_nearest_transcript_end(algn_pair[0].iv.end_d, gene_ids_union, transcript_ends, is_randomn)
        counter[gene_id] += 1
        return
    
    # both intersect and union are empty
    # Gene Match (introns)
    gene_ids_union, gene_ids_intersect = get_union_intersect_genes(algn_pair, genes)
    # If still not found, try flanking regions (after polyA cleavage point)
    if len(gene_ids_intersect) == 0 and len(gene_ids_union) == 0:
        gene_ids_union, gene_ids_intersect = get_union_intersect_flanking_genes(algn_pair, genes)
    # If still not found, try the opposite strand
    if len(gene_ids_intersect) == 0 and len(gene_ids_union) == 0:
        reverse_strand(algn_pair[0])
        reverse_strand(algn_pair[1])
        gene_ids_union, gene_ids_intersect = get_union_intersect_genes(algn_pair, exons)
        if len(gene_ids_intersect) > 0:
            suffix = "" # remove _intron suffix (because we found on the opposite strand by exons)
    if len(gene_ids_intersect) == 0 and len(gene_ids_union) == 0:
        gene_ids_union, gene_ids_intersect = get_union_intersect_genes(algn_pair, genes)
    
    if len(gene_ids_intersect) == 1:
        counter[next(iter(gene_ids_intersect)) + suffix] += 1
    elif len(gene_ids_intersect) > 1:
        gene_id = find_nearest_transcript_end(algn_pair[0].iv.end_d, gene_ids_intersect, transcript_ends, is_randomn)
        counter[gene_id + suffix] += 1
    elif len(gene_ids_union) == 1:
        counter[next(iter(gene_ids_union)) + suffix] += 1
    elif len(gene_ids_union) > 1:
        gene_id = find_nearest_transcript_end(algn_pair[0].iv.end_d, gene_ids_union, transcript_ends, is_randomn)
        counter[gene_id + suffix] += 1
    else: # No gene found
        counter["_no_feature"] += 1

def get_union_intersect_genes(algn_pair, htseq_array)-> tuple[set, set]:
    """
    Get the union and intersection of genes from a pair of alignments

    htseq_array: HTSeq.GenomicArrayOfSets (genes or exons)
    """
    genes_union = set()
    genes_intersect = set()
    is_intersect_empty = True
    for op in algn_pair[0].cigar:
        if op.type != "M":
            continue
        for iv, gene in htseq_array[op.ref_iv].steps():
            genes_union |= gene
            if is_intersect_empty:
                genes_intersect |= gene
                is_intersect_empty = False
            else:
                genes_intersect &= gene
    
    for op in algn_pair[1].cigar:
        if op.type != "M":
            continue
        for iv, gene in htseq_array[op.ref_iv].steps():
            genes_union |= gene
            if is_intersect_empty:
                genes_intersect |= gene
                is_intersect_empty = False
            else:
                genes_intersect &= gene
    return genes_union, genes_intersect

def find_nearest_transcript_end(r1_end, gene_ids_intersect, transcript_ends, is_randomn, region_size=100)-> str:
    """
    Find the gene name of nearest transcript end to the given position (short-dT only)

    And make sure there is no another gene within ± 100 bp of the nearest transcript end

    transcript_ends: dict { gene_name: numpy.array of transcript ends }
    """
    gene_id = "_ambiguous"

    if is_randomn:
        return gene_id # _ambiguous for randomN RT barcodes

    gene_id_ends = dict()
    for gene_id in gene_ids_intersect:
        if gene_id in transcript_ends:
            gene_id_ends[gene_id] = abs(transcript_ends[gene_id] - r1_end).min()
        else:
            print("WARNING: Gene not found in transcript ends", gene_id)

    gene_end_min = min(gene_id_ends.values())

    is_more_than_one = False
    for gene_id in gene_id_ends:
        if gene_id_ends[gene_id] < gene_end_min + region_size:
            if is_more_than_one: # More than one gene within 100 bp
                break
            else:
                gene_id = gene_id
                is_more_than_one = True
    
    return gene_id

def get_union_intersect_flanking_genes(algn_pair, htseq_array, flanking_size=1000)-> tuple[set, set]:
    """
    Get genes count from 3' flanking regions using a pair of alignments
    (The reads may come from pre-mRNA after cleavage point AAUAAA)

    On + strand => interval - flanking_size (need to make sure start is not negative)
    On - strand => interval + flanking_size
    """
    gene_ids_union = set()
    gene_ids_intersect = set()
    is_intersect_empty = True

    if algn_pair[1].iv.strand == "+":
        flanking_size = -flanking_size # for + strand, we need to subtract

    diff = flanking_size

    # make sure start is not negative
    if diff < 0 and algn_pair[1].iv.start <= -diff:
        diff = -algn_pair[1].iv.start

    for op in algn_pair[1].cigar:
        if op.type != "M":
            continue
        
        op.ref_iv.start += diff
        op.ref_iv.end += diff

        for iv, gene in htseq_array[op.ref_iv].steps():
            gene_ids_union |= gene
            if is_intersect_empty:
                gene_ids_intersect |= gene
                is_intersect_empty = False
            else:
                gene_ids_intersect &= gene
        
        op.ref_iv.start -= diff
        op.ref_iv.end -= diff
    
    if diff < 0 and algn_pair[0].iv.start <= -diff:
        diff = -algn_pair[0].iv.start
    
    for op in algn_pair[0].cigar:
        if op.type != "M":
            continue

        op.ref_iv.start += diff
        op.ref_iv.end += diff

        for iv, gene in htseq_array[op.ref_iv].steps():
            gene_ids_union |= gene
            if is_intersect_empty:
                gene_ids_intersect |= gene
                is_intersect_empty = False
            else:
                gene_ids_intersect &= gene

        op.ref_iv.start -= diff
        op.ref_iv.end -= diff

    return gene_ids_union, gene_ids_intersect

def reverse_strand(algn):
    """
    Reverse the strand of an alignment and its cigar operations' reference intervals (in-place)
    """
    # algn is reference, so we can change it in-place, no need to return
    algn.iv.strand = "-" if algn.iv.strand == "+" else "+"
    for op in algn.cigar:
        op.ref_iv.strand = algn.iv.strand

def count_bam(bam_path, output_prefix, genes, exons, transcript_ends, gene_id_index, randomn_barcodes):
    """
    Count gene reads in a name sorted BAM file
    """
    alignments = HTSeq.BAM_Reader(bam_path)
    output_count_mtx = gzip.open(f"{output_prefix}.count.gz", "wt") if IS_GZIP else open(f"{output_prefix}.count", "w")
    output_count_mtx.write("gene_index,cell_index,count\n")
    output_cell_ids = gzip.open(f"{output_prefix}_cell_ids.csv.gz", "wt") if IS_GZIP else open(f"{output_prefix}_cell_ids.csv", "w")
    output_cell_ids.write("cell_index,cell_id,unmatched_rate\n")
    i7_prefix = os.path.basename(output_prefix)

    previous_barcode = ""
    cell_index = 0 # cell index starts from 0 for later sparse matrix creation
    counter = collections.Counter()

    for bundle in HTSeq.pair_SAM_alignments(alignments, bundle=True):
        if len(bundle) != 1:
            continue

        algn_pair = bundle[0]
        # Since the bam file is name sorted and the name starts with cell barcode, we can use it to separate cells
        # NAME:@<cell_barcode>,<UMI>,<read_name>
        barcode = algn_pair[0].read.name.partition(",")[0] # @<ligation_barcode>-<RT_barcode>,UMI,original_read_name
        barcode_ligation, barcode_rt = barcode.split("-") # @<ligation_barcode>-<RT_barcode>
        is_randomn = barcode_rt in randomn_barcodes

        if barcode == previous_barcode: # Same cell
            count_pair(algn_pair, genes, exons, transcript_ends, is_randomn, counter)
        else: # New cell
            if len(counter) == 0:
                count_pair(algn_pair, genes, exons, transcript_ends, is_randomn, counter)
                previous_barcode = barcode
            else: # Save the counts and reset
                # record stats of current cell
                record_results(counter, output_count_mtx, output_cell_ids, 
                               cell_index, i7_prefix, previous_barcode, gene_id_index)
                # reset
                counter = collections.Counter()
                cell_index += 1
                count_pair(algn_pair, genes, exons, transcript_ends, is_randomn, counter)
                previous_barcode = barcode
    
    # record the last cell
    record_results(counter, output_count_mtx, output_cell_ids, 
                   cell_index, i7_prefix, previous_barcode, gene_id_index)
    
    alignments.close()
    output_count_mtx.close()
    output_cell_ids.close()

def record_results(counter, output_count_mtx, output_cell_ids, cell_index, i7_prefix, barcode, gene_id_index):
    """
    Record the results of the current cell
    """
    n_reads = n_unmatched = n_ambiguous = 0
    for exon in counter:
        if exon.startswith("_"): # {"_unmapped", "_ambiguous", "_ambiguous_intron", "_no_feature"}
            if exon == "_no_feature":
                n_unmatched += counter[exon]
            else:
                n_ambiguous += counter[exon]
        else:
            n_reads += counter[exon]
            # gene_index, cell_index, count (index starts from 0 for later sparse matrix creation)
            output_count_mtx.write(f"{gene_id_index[exon]},{cell_index},{counter[exon]}\n")
    
    unmatched_rate = n_unmatched / (n_unmatched + n_reads + n_ambiguous)
    ligation_barcode, rt_barcode = barcode.split("-")
    output_cell_ids.write(f"{cell_index},{i7_prefix}.{ligation_barcode}.{rt_barcode},{unmatched_rate}\n") # <i7_prefix>-<ligation_barcode>-<RT_barcode>

def count_bam_parallel(input_folder, output_folder, genes, exons, transcript_ends, gene_id_index, randomn_barcodes, i7_prefix):
    """
    Count exon reads in BAM files for a i7 file prefix (used for parallel processing)
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {i7_prefix} Running...")
    bam_path = os.path.join(input_folder, f"{i7_prefix}.filtered.dedup.bam")
    output_prefix = os.path.join(output_folder, i7_prefix)
    count_bam(bam_path, output_prefix, genes, exons, transcript_ends, gene_id_index, randomn_barcodes)
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {i7_prefix} Done")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Count exon expression per cell")
    parser.add_argument("-i", "--input_folder", help="Input folder with BAM files", required=True)
    parser.add_argument("-o", "--output_folder", help="Output folder for count files", required=True)
    parser.add_argument("-p", "--i7_prefix_file", help="Text file with i7 demultiplexed prefix", required=True)
    parser.add_argument("-g", "--gtf", help="GTF file with exon information", required=True)
    parser.add_argument("-rt", "--rt_barcode_tsv", help="RT barcode tsv (Sample	shortdt	randomN)", required=True)
    parser.add_argument("--gzip", action="store_true", help="Compress the output files")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads")
    args = parser.parse_args()
    
    file_check(args.gtf) # safety check
    file_check(args.i7_prefix_file)
    file_check(args.rt_barcode_tsv)
    os.makedirs(args.output_folder, exist_ok=True)
    global IS_GZIP
    IS_GZIP = args.gzip

    genes, exons, gene_transcript_ends, genes_table = get_genes_exons_array(args.gtf)
    out_gene_annotation_path = os.path.join(args.output_folder, "gene_ids.csv")
    if IS_GZIP:
        out_gene_annotation_path += ".gz"
    genes_table.to_csv(out_gene_annotation_path, index=False)
    
    gene_id_index = genes_table.set_index("gene_id")["gene_index"].to_dict()
    del genes_table # free memory
    randomn_barcodes = set(pd.read_csv(args.rt_barcode_tsv, sep="\t")["randomN"])
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] GTF and RT barcodes loaded")

    from functools import partial
    from multiprocessing import Pool
    i7_prefixes = list(prefix for prefix in open(args.i7_prefix_file).read().splitlines() if prefix)
    worker_count_bam = partial(
        count_bam_parallel, args.input_folder, args.output_folder,
        genes, exons, gene_transcript_ends, gene_id_index, randomn_barcodes)
    with Pool(args.threads) as pool:
        pool.map(worker_count_bam, i7_prefixes)

def file_check(file_path: str):
    """
    Check if the file exists, otherwise exit(1)
    """
    if not os.path.exists(file_path):
        print(f"Error: {file_path} does not exist")
        exit(1)

def get_genes_exons_array(gtf_path: str):
    """
    Get the genes and exons array from a GTF file
    Also get the gene_transcript_ends and gene_count_table
    """
    genes = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    gene_transcript_ends = {}

    gene_count_index = 0 # starts from 0 for later sparse matrix creation
    indexes = []
    gene_ids = []
    gene_types = []
    gene_names = []
    exon_introns = []

    for feature in HTSeq.GFF_Reader(gtf_path, end_included=True):
        if feature.type == "exon":
            exons[feature.iv] += feature.attr["gene_id"]
        elif feature.type == "transcript":
            if feature.attr["gene_id"] in gene_transcript_ends:
                gene_transcript_ends[feature.attr["gene_id"]].add(feature.iv.end_d)
            else:
                gene_transcript_ends[feature.attr["gene_id"]] = set([feature.iv.end_d])
        elif feature.type == "gene":
            genes[feature.iv] += feature.attr["gene_id"]
            indexes.append(gene_count_index)
            gene_ids.append(feature.attr["gene_id"])
            gene_types.append(feature.attr["gene_type"])
            gene_names.append(feature.attr["gene_name"])
            exon_introns.append("exon")
            gene_count_index += 1

            # intron counting table
            indexes.append(gene_count_index)
            gene_ids.append(f"{feature.attr['gene_id']}_intron")
            gene_types.append(feature.attr["gene_type"])
            gene_names.append(f"{feature.attr['gene_name']}_intron")
            exon_introns.append("intron")
            gene_count_index += 1

    gene_count_table = pd.DataFrame({
        "gene_index": indexes,
        "gene_id": gene_ids,
        "gene_type": gene_types,
        "gene_name": gene_names,
        "feature": exon_introns
    })

    # Convert set to numpy array for faster ends checking
    for gene_id in gene_transcript_ends:
        gene_transcript_ends[gene_id] = np.array(list(gene_transcript_ends[gene_id]))

    return genes, exons, gene_transcript_ends, gene_count_table

if __name__ == "__main__":
    main()
