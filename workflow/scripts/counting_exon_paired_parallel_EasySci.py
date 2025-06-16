#!/usr/bin/env python3

##############################################################
# counting_exon_paired_parallel.py
# Description: This script attaches ligation and RT barcodes and UMI to FASTQ files
#   Input:
#      --input_folder:   Folder with BAM files
#      --output_folder:  Output folder for count files
#      --i7_prefix_file: Text file with i7 demultiplexed prefix
#      --gtf:            GTF file with exon information only
#      --threads:        Number of threads
#      --gzip:           Compress the output files
#   Output:
#      <output_folder>/<i7_prefix>.count.gz: Count file with gene expression per cell ┐ 
#      <output_folder>/<i7_prefix>_cell_id.csv.gz: Cell id file with unmatched rate   │ for later sparse matrix creation
#      <output_folder>/gene_exon_ids.csv.gz: Gene exon id table                       ┘
##############################################################

# Counting Strategy:
# 1. Try to match two reads to exons of the same gene (stranded)
# 2. If not found on expected strand, try the opposite strand
# 3. ambigous if multiple genes are found
# 4. When multiple exons from the same gene are found, number will be 0.5/n_exons for each exon

# If only one of the reads is mapped to exons, while the other is not, then:
#     If mapped read is uniquely mapped to a gene
#         True -> count exons
#         False -> ambiguous
# If both reads are mapped to only 1 exon, then:
#     If the exons are from the same gene
#         True -> count both
#         False -> ambiguous
# If both reads are mapped to multiple exons, then:
#     Try to find the intersection gene (only one)
#         True -> count exons from the intersection gene
#         False -> ambiguous
# If one read is mapped to multiple exons, while the other is mapped to only 1 exon, then:
#     For the only 1 exon read, count it
#     For the multiple exons read, only count exons from the same gene

import collections
from datetime import datetime
import HTSeq
import os
import pandas as pd
import gzip
IS_GZIP = False

def count_pair(algn_pair, exons, counter):
    """
    Count exon reads in a pair of alignments
    """
    if not algn_pair[0].aligned or not algn_pair[1].aligned:
        counter["_unmapped"] += 1 # Just in case, since we only have primary alignments
        return
        # Because Counter is a reference here, we don't need to return it.
    if algn_pair[0].iv.chrom not in exons.chrom_vectors or algn_pair[1].iv.chrom not in exons.chrom_vectors:
        counter["_unmapped"] += 1 # Just in case, since we use the same GTF file for alignment and counting
        return

    # Because of the Illumina sequencing and library preparation.
    #   1st read is always anti-sense to the gene(mRNA), while 2nd read is sense.
    reverse_strand(algn_pair[0])
    first_algn_exons = get_exon_reads_exons(algn_pair[0], exons)
    second_algn_exons = get_exon_reads_exons(algn_pair[1], exons)

    # But also because this is single-nucleus RNA-seq, must have pre-mRNA.
    #   So, it is worth checking the opposite strand mapping.
    if len(first_algn_exons) == 0 and len(second_algn_exons) == 0:
        reverse_strand(algn_pair[0])
        first_algn_exons = get_exon_reads_exons(algn_pair[0], exons)
        reverse_strand(algn_pair[1])
        second_algn_exons = get_exon_reads_exons(algn_pair[1], exons)

    # If still not matched, then _no_feature
    if len(first_algn_exons) == 0 and len(second_algn_exons) == 0:
        counter["_no_feature"] += 1

    elif len(first_algn_exons) == 0: # len(second_algn_exons) > 0
        count_yes_no(counter, second_algn_exons)

    elif len(second_algn_exons) == 0: # len(first_algn_exons) > 0
        count_yes_no(counter, first_algn_exons)

    elif len(first_algn_exons) == 1 and len(second_algn_exons) == 1:
        count_one_one(counter, next(iter(first_algn_exons)), next(iter(second_algn_exons)))
                      
    elif len(first_algn_exons) > 1 and len(second_algn_exons) > 1:
        count_multiple_multiple(counter, first_algn_exons, second_algn_exons)

    elif len(first_algn_exons) > 1: # len(second_algn_exons) == 1
        count_one_multiple(counter, next(iter(second_algn_exons)), first_algn_exons)

    else: # len(first_algn_exons) == 1 and len(second_algn_exons) > 1
        count_one_multiple(counter, next(iter(first_algn_exons)), second_algn_exons)

def count_yes_no(counter, overlap_exons: set):
    """
    Use this function when one read has overlapping exons while the other has no overlap

    - overlap_exons: set of gene_exon_id that the matched alignment is mapped to

    If mapped read is uniquely mapped to a gene, count exons (1/n_exons for each exon)
    """
    counter["_no_feature"] += 0.5
    algn_genes = set(exon.partition("-")[0] for exon in overlap_exons)
    if len(algn_genes) == 1:
        for exon in overlap_exons:
            counter[exon] += 0.5 / len(overlap_exons)
    else:
        counter["_ambiguous"] += 0.5

def count_one_one(counter, overlap_exon1: str, overlap_exon2: str):
    """
    Use this function when both reads have only 1 overlapping exon

    - overlap_exon1: gene_exon_id that the 1st read is mapped to
    - overlap_exon2: gene_exon_id that the 2nd read is mapped to

    If both reads are mapped to the same gene, then count them separately

    next(iter(set)) to get the only element of a set
    """
    if overlap_exon1.partition("-")[0] == overlap_exon2.partition("-")[0]:
        counter[overlap_exon1] += 0.5
        counter[overlap_exon2] += 0.5
    else:
        counter["_ambiguous"] += 1

def count_multiple_multiple(counter, overlap_exons1: set, overlap_exons2: set):
    """
    Use this function when both reads have multiple overlapping exons

    - overlap_exons1: set of gene_exon_id that the 1st read is mapped to
    - overlap_exons2: set of gene_exon_id that the 2nd read is mapped to

    If both reads have uniquely intersected genes, then count them separately
        only count exons belonging to the intersection gene (1/n_exons for each exon)
    """
    algn_genes1 = set(exon.partition("-")[0] for exon in overlap_exons1)
    algn_genes2 = set(exon.partition("-")[0] for exon in overlap_exons2)
    intersection_genes = algn_genes1 & algn_genes2
    if len(intersection_genes) == 1:
        intersection_gene = next(iter(intersection_genes))
        overlap_exons1 = [exon for exon in overlap_exons1 if exon.startswith(intersection_gene)]
        for exon in overlap_exons1:
            counter[exon] += 0.5 / len(overlap_exons1)
        overlap_exons2 = [exon for exon in overlap_exons2 if exon.startswith(intersection_gene)]
        for exon in overlap_exons2:
            counter[exon] += 0.5 / len(overlap_exons2)
    else:
        counter["_ambiguous"] += 1

def count_one_multiple(counter, overlap_exon: str, overlap_exons: set):
    """
    Use this function when one read has only 1 overlapping exon while the other has multiple overlapping exons

    - overlap_exon: gene_exon_id that the 1st read is mapped
    - overlap_exons: set of gene_exon_id that the 2nd read is mapped to

    Counting
    - exon: overlap_exon is counted ???
    - set:  only count exons belonging to the same gene as overlap_exon (1/n_exons for each exon)

    next(iter(set)) to get the only element of a set
    """
    overlap_gene = overlap_exon.partition("-")[0]
    overlap_genes = set(exon.partition("-")[0] for exon in overlap_exons)

    counter[overlap_exon] += 0.5
    if overlap_gene in overlap_genes:
        overlap_exons = [exon for exon in overlap_exons if exon.startswith(overlap_gene)]
        for exon in overlap_exons:
            counter[exon] += 0.5 / len(overlap_exons)
    else:
        counter["_ambiguous"] += 0.5

def get_exon_reads_exons(algn, exons)-> set:
    """
    Only count reads in exons (no introns or junctions)
    """
    algn_exons = set()
    first_run = True
    for op in algn.cigar:
        if op.type != "M":
            continue
        for iv, exon in exons[op.ref_iv].steps():
            if first_run:
                algn_exons |= exon
                first_run = False
            else:
                algn_exons &= exon
    return algn_exons

def reverse_strand(algn):
    """
    Reverse the strand of an alignment and its cigar operations' reference intervals (in-place)
    """
    # algn is reference, so we can change it in-place, no need to return
    algn.iv.strand = "-" if algn.iv.strand == "+" else "+"
    for op in algn.cigar:
        op.ref_iv.strand = algn.iv.strand

def count_bam(bam_path, output_prefix, exons, gene_exon_id_index):
    """
    Count exon reads in a name sorted BAM file

    exons: HTSeq.GenomicArrayOfSets of all gene_id-exon_id
    gene_exon_id_index: { gene_exon_id: exon_index }
    """
    alignments = HTSeq.BAM_Reader(bam_path)
    output_count_mtx = gzip.open(f"{output_prefix}.count.gz", "wt") if IS_GZIP else open(f"{output_prefix}.count", "w")
    output_count_mtx.write("exon_index,cell_index,count\n")
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
        barcode = algn_pair[0].read.name.partition(",")[0]

        if barcode == previous_barcode: # Same cell
            count_pair(algn_pair, exons, counter)
        else: # New cell
            if len(counter) == 0:
                count_pair(algn_pair, exons, counter)
                previous_barcode = barcode
            else: # Save the counts and reset
                # record stats of current cell
                record_results(counter, output_count_mtx, output_cell_ids, 
                               cell_index, i7_prefix, previous_barcode, gene_exon_id_index)
                # reset
                counter = collections.Counter()
                cell_index += 1
                count_pair(algn_pair, exons, counter)
                previous_barcode = barcode
    
    # record the last cell
    record_results(counter, output_count_mtx, output_cell_ids, 
                   cell_index, i7_prefix, previous_barcode, gene_exon_id_index)
    
    alignments.close()
    output_count_mtx.close()
    output_cell_ids.close()

def record_results(counter, output_count_mtx, output_cell_ids, cell_index, i7_prefix, barcode, gene_exon_id_index):
    """
    Record the results of the current cell
    """
    n_reads = n_unmatched = n_ambiguous = 0
    for exon in counter:
        if exon.startswith("_"): # {"_unmapped", "_no_feature", "_ambiguous"}
            if exon == "_no_feature":
                n_unmatched += counter[exon]
            else:
                n_ambiguous += counter[exon]
        else:
            n_reads += counter[exon]
            # exon_index, cell_index, count (index starts from 0 for later sparse matrix creation)
            output_count_mtx.write(f"{gene_exon_id_index[exon]},{cell_index},{counter[exon]}\n")
    
    unmatched_rate = n_unmatched / (n_unmatched + n_reads + n_ambiguous)
    output_cell_ids.write(f"{cell_index},{i7_prefix}.{barcode},{unmatched_rate}\n")

def count_bam_parallel(input_folder, output_folder, exons, gene_exon_id_index, i7_prefix):
    """
    Count exon reads in BAM files for a i7 file prefix (used for parallel processing)
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {i7_prefix} Running...")
    bam_path = os.path.join(input_folder, f"{i7_prefix}.filtered.dedup.bam")
    output_prefix = os.path.join(output_folder, i7_prefix)
    count_bam(bam_path, output_prefix, exons, gene_exon_id_index)
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {i7_prefix} Done")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Count exon expression per cell")
    parser.add_argument("-i", "--input_folder", help="Input folder with BAM files", required=True)
    parser.add_argument("-o", "--output_folder", help="Output folder for count files", required=True)
    parser.add_argument("-p", "--i7_prefix_file", help="Text file with i7 demultiplexed prefix", required=True)
    parser.add_argument("-g", "--gtf", help="GTF file with exon information", required=True)
    parser.add_argument("--gzip", action="store_true", help="Compress the output files")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads")
    args = parser.parse_args()

    file_check(args.gtf) # safety check
    file_check(args.i7_prefix_file)
    os.makedirs(args.output_folder, exist_ok=True)
    global IS_GZIP
    IS_GZIP = args.gzip

    exons, exons_table = get_exons_array(args.gtf)
    out_exon_annotation_path = os.path.join(args.output_folder, "gene_exon_ids.csv.gz" if IS_GZIP else "gene_exon_ids.csv")
    exons_table.to_csv(out_exon_annotation_path, index=False)

    gene_exon_id_index = exons_table.set_index("gene_exon_id")["exon_index"].to_dict()
    del exons_table # free memory

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Exon reference loaded")
    from functools import partial
    from multiprocessing import Pool
    i7_prefixes = list(prefix for prefix in open(args.i7_prefix_file).read().splitlines() if prefix)
    worker_count_bam = partial(count_bam_parallel, args.input_folder, args.output_folder, exons, gene_exon_id_index)
    with Pool(args.threads) as pool:
        pool.map(worker_count_bam, i7_prefixes)

def file_check(file_path: str):
    """
    Check if the file exists, otherwise exit(1)
    """
    if not os.path.exists(file_path):
        print(f"Error: {file_path} does not exist")
        exit(1)

def get_exons_array(gtf_path: str) -> tuple[HTSeq.GenomicArrayOfSets, pd.DataFrame]:
    """
    Get the exons array and gene_exon_id table from the GTF file
    """
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    exon_count_index = 0 # starts from 0 for later sparse matrix creation
    gene_exon_start_end_set = set() # to check if the exon is duplicated
                                    # format (gene_id-start-end) since gene can only on one strand of chromosome
                                    # If the same gene has multiple exons with the same start and end, only keep the first exon
    gene_exon_id_set = set() # set check is faster than list
    indexes = []
    gene_exon_ids = []
    gene_types = []
    gene_names = []
    chroms = []
    strands = []
    starts = []
    ends = []
    lengths = []

    for feature in HTSeq.GFF_Reader(gtf_path):
        if feature.type == "exon":
            gene_exon_id = f"{feature.attr['gene_id']}-{feature.attr['exon_id']}"
            if gene_exon_id in gene_exon_id_set: # Skip duplicated exons 
                # (different trnascripts from the same gene may have the same exon)
                continue
            gene_start_end = f"{feature.attr['gene_id']}-{feature.iv.start}-{feature.iv.end}"
            if gene_start_end in gene_exon_start_end_set: # Skip duplicated exons
                # (different trnascripts from the same gene may have the same exon)
                continue
            gene_exon_start_end_set.add(gene_start_end)

            exons[feature.iv] += gene_exon_id
            gene_exon_id_set.add(gene_exon_id)

            indexes.append(exon_count_index)
            gene_exon_ids.append(gene_exon_id)
            gene_types.append(feature.attr['gene_type'])
            gene_names.append(feature.attr['gene_name'])
            chroms.append(feature.iv.chrom)
            strands.append(feature.iv.strand)
            lengths.append(feature.iv.end - feature.iv.start) # Length of the exon (be careful with 0-based indexing)
            starts.append(feature.iv.start + 1)  # HTSeq uses 0-based indexing, but we use 1-based for output
            ends.append(feature.iv.end) # HTSeq uses 0-based indexing but the end is exclusive, so we can use it directly.
            exon_count_index += 1

    exon_count_table = pd.DataFrame({
        "exon_index": indexes,
        "gene_exon_id": gene_exon_ids,
        "gene_type": gene_types,
        "gene_name": gene_names,
        "chrom": chroms,
        "strand": strands,
        "start": starts,
        "end": ends,
        "length": lengths
    })

    return exons, exon_count_table

if __name__ == "__main__":
    main()
