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
# Make sure PE reads are in the same gene. If multiple genes are found, then ambiguous.
# Exon reads have only one block.
# Junction reads have multiple blocks, which is separated by "N" in the CIGAR string.

import collections
from datetime import datetime
import HTSeq
import os
import pandas as pd
import gzip
IS_GZIP = False

def count_pair_with_junctions(algn_pair, exons, exon_gene_dict, counter):
    """
    Count exon reads in a pair of alignments with junctions
    """
    if not algn_pair[0].aligned or not algn_pair[1].aligned:
        counter["_unmapped"] += 1 # Just in case, since we only have primary alignments
        return
        # Because Counter is a reference here, we don't need to return it.

    if algn_pair[0].iv.chrom not in exons.chrom_vectors or algn_pair[1].iv.chrom not in exons.chrom_vectors:
        counter["_unmapped"] += 1 # Just in case, since we use the same GTF file for alignment and counting
        return

    reverse_strand(algn_pair[0]) # You need to do this because of the Illumina sequencing and library preparation.

    # Make sure first alignment is the one that starts first on the genome
    if algn_pair[0].iv.start < algn_pair[1].iv.start:
        first_algn, second_algn = algn_pair[0], algn_pair[1]
    else:
        first_algn, second_algn = algn_pair[1], algn_pair[0]

    first_algn_blocks_and_genes = get_exon_junction_reads_blocks(first_algn, exons, exon_gene_dict)
    second_algn_blocks_and_genes = get_exon_junction_reads_blocks(second_algn, exons, exon_gene_dict)
    # The blocks is a tuple of:
    # - [0] exon blocks: list of sets where each set contains gene_exon_id from the same block
    #   - Each block is separated by "N" in the CIGAR string.
    # - [1] intersection genes: set of genes that all blocks have in common (skip empty blocks)
    # - [2] union genes: set of all genes that are present in any of the blocks

    # any to all
    if all(len(block) == 0 for block in first_algn_blocks_and_genes[0]) and all(len(block) == 0 for block in second_algn_blocks_and_genes[0]):
        # No feature found in both reads, try the opposite strand
        reverse_strand(first_algn)
        first_algn_blocks_and_genes = get_exon_junction_reads_blocks(first_algn, exons, exon_gene_dict)
        reverse_strand(second_algn)
        second_algn_blocks_and_genes = get_exon_junction_reads_blocks(second_algn, exons, exon_gene_dict)

    # If still not matched, then _no_feature
    if all(len(block) == 0 for block in first_algn_blocks_and_genes[0]) and all(len(block) == 0 for block in second_algn_blocks_and_genes[0]):
        counter["_no_feature"] += 1
        return

    elif len(first_algn_blocks_and_genes[0]) == 0: # len(second_algn_blocks_and_genes[0]) > 0
        if len(second_algn_blocks_and_genes[1]) == 1: # only one gene
            count_exon_blocks_using_gene(counter, exon_gene_dict, next(iter(second_algn_blocks_and_genes[1])), second_algn_blocks_and_genes[0])
        elif len(second_algn_blocks_and_genes[2]) == 1:
            count_exon_blocks_using_gene(counter, exon_gene_dict, next(iter(second_algn_blocks_and_genes[2])), second_algn_blocks_and_genes[0])
        else:
            counter["_ambiguous"] += 1
        return

    elif len(second_algn_blocks_and_genes[0]) == 0: # len(first_algn_blocks_and_genes[0]) > 0
        if len(first_algn_blocks_and_genes[1]) == 1: # only one gene
            count_exon_blocks_using_gene(counter, exon_gene_dict, next(iter(first_algn_blocks_and_genes[1])), first_algn_blocks_and_genes[0])
        elif len(first_algn_blocks_and_genes[2]) == 1:
            count_exon_blocks_using_gene(counter, exon_gene_dict, next(iter(first_algn_blocks_and_genes[2])), first_algn_blocks_and_genes[0])
        else:
            counter["_ambiguous"] += 1
        return
    
    # If both reads have more than 0 blocks, then we need to check the intersection of genes
    common_genes = first_algn_blocks_and_genes[1] & second_algn_blocks_and_genes[1]
    if len(common_genes) != 1:
        common_genes = first_algn_blocks_and_genes[1] & second_algn_blocks_and_genes[2]
    if len(common_genes) != 1:
        common_genes = first_algn_blocks_and_genes[2] & second_algn_blocks_and_genes[1]
    if len(common_genes) != 1:
        common_genes = first_algn_blocks_and_genes[2] & second_algn_blocks_and_genes[2]
    if len(common_genes) != 1:
        counter["_ambiguous"] += 1
        return
    
    # Try to see there are overlap blocks
    # If first end > second start, then they overlap
    blocks_to_count = []
    if first_algn.iv.end <= second_algn.iv.start:
        # No overlap, try to see if there is any common exons on the edges
        common_exons = first_algn_blocks_and_genes[0][-1] & second_algn_blocks_and_genes[0][0]
        if common_exons:
            blocks_to_count.append(common_exons)
            # Skip the last block of the first read and the first block of the second read
            if len(first_algn_blocks_and_genes[0]) > 1:
                blocks_to_count.extend(first_algn_blocks_and_genes[0][:-1])
            if len(second_algn_blocks_and_genes[0]) > 1:
                blocks_to_count.extend(second_algn_blocks_and_genes[0][1:])
        else:
            blocks_to_count.extend(first_algn_blocks_and_genes[0])
            blocks_to_count.extend(second_algn_blocks_and_genes[0])
    else:
        skip_indices_first = []
        skip_indices_second = []
        
        # # Try to identify intersecting blocks
        # i = len(first_algn_blocks_and_genes[0]) - 1 # last block of the first read
        # j = len(second_algn_blocks_and_genes[0]) - 1 # last block of the second read
        # is_common_found = False
        # while i >= 0 and j >= 0:
        #     common_exons = first_algn_blocks_and_genes[0][i] & second_algn_blocks_and_genes[0][j]
        #     if common_exons:
        #         blocks_to_count.append(common_exons)
        #         is_common_found = True
        #         # Skip the block of the first read and the second read
        #         skip_indices_first.append(i)
        #         skip_indices_second.append(j)
        #     if is_common_found:
        #         i -= 1
        #         j -= 1
        #     else:
        #         # If no common exons, then check the next previous block of the second read
        #         j -= 1

        for i in range(len(first_algn_blocks_and_genes[0])):
            for j in range(len(second_algn_blocks_and_genes[0])):
                common_exons = first_algn_blocks_and_genes[0][i] & second_algn_blocks_and_genes[0][j]
                if common_exons:
                    blocks_to_count.append(common_exons)
                    # Skip the block of the first read and the second read
                    skip_indices_first.append(i)
                    skip_indices_second.append(j)
        
        blocks_to_count.extend(
            first_algn_blocks_and_genes[0][i] for i in range(len(first_algn_blocks_and_genes[0])) if i not in skip_indices_first
        )
        blocks_to_count.extend(
            second_algn_blocks_and_genes[0][j] for j in range(len(second_algn_blocks_and_genes[0])) if j not in skip_indices_second
        )

    count_exon_blocks_using_gene(counter, exon_gene_dict, next(iter(common_genes)), blocks_to_count)

def get_exon_junction_reads_blocks(algn, exons, exon_gene_dict)-> tuple[list[set], set, set]:
    """
    For a given alignment, return a list of exon sets.
    Each set corresponds to a continuous aligned block (separated by N in the CIGAR string).

    - algn: HTSeq.Alignment object
    - exons: HTSeq.GenomicArrayOfSets of all gene_id-exon_id

    return a tuple of:
    - [0] exon blocks: list of sets where each set contains gene_exon_id from the same block
        - Each block is separated by "N" in the CIGAR string.
    - [1] intersection genes: set of genes that all blocks have in common (skip empty blocks)
    - [2] union genes: set of all genes that are present in any of the blocks
    """
    exon_blocks = []
    intersection_genes = set()
    union_genes = set()

    current_exons = set()
    is_first_op = True

    for op in algn.cigar:
        if op.type == "N":  # Spliced region (e.g., intron)
            exon_blocks.append(current_exons)
            # Reset for the next block
            current_exons = set()
            is_first_op = True
            continue

        if op.type != "M":
            continue

        for iv, exon_ids in exons[op.ref_iv].steps():
            if is_first_op:
                current_exons |= exon_ids
                is_first_op = False
            else:
                current_exons &= exon_ids

    exon_blocks.append(current_exons)

    # Get genes from the blocks (except empty blocks)
    blocks_genes = [set(exon_gene_dict[exon] for exon in block) for block in exon_blocks if block]
    if blocks_genes:
        intersection_genes = set.intersection(*blocks_genes)
        union_genes = set.union(*blocks_genes)

    return (exon_blocks, intersection_genes, union_genes)

def count_exon_blocks_using_gene(counter, exon_gene_dict, gene, exon_blocks: list[set], multiplier = 1):
    """
    Count exon blocks using the gene name
    - counter: collections.Counter object to store counts
    - exon_gene_dict: dict mapping gene_exon_id to gene_id
    - gene: gene name (e.g., "ENSG00000123456")
    - exon_blocks: list of sets where each set contains gene_exon_id from the same block
    """
    is_counted = False

    for block in exon_blocks:
        block_exons_from_the_same_gene = [exon for exon in block if exon_gene_dict[exon] == gene]
        if block_exons_from_the_same_gene:
            for exon in block_exons_from_the_same_gene:
                is_counted = True
                counter[exon] += multiplier / len(block_exons_from_the_same_gene)
    
    if not is_counted:
        counter["_ambiguous"] += 1  # If no exons from the same gene are found, count as ambiguous

def reverse_strand(algn):
    """
    Reverse the strand of an alignment and its cigar operations' reference intervals (in-place)
    """
    # algn is reference, so we can change it in-place, no need to return
    algn.iv.strand = "-" if algn.iv.strand == "+" else "+"
    for op in algn.cigar:
        op.ref_iv.strand = algn.iv.strand

def count_bam(bam_path, output_prefix, exons, gene_exon_id_index, exon_gene_dict):
    """
    Count exon reads in a name sorted BAM file

    exons: HTSeq.GenomicArrayOfSets of all gene_id-exon_id
    gene_exon_id_index: { gene_exon_id: exon_index }
    exon_gene_dict: { gene_exon_id: gene_id }
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
            count_pair_with_junctions(algn_pair, exons, exon_gene_dict, counter)
        else: # New cell
            if len(counter) == 0:
                count_pair_with_junctions(algn_pair, exons, exon_gene_dict, counter)
                previous_barcode = barcode
            else: # Save the counts and reset
                # record stats of current cell
                record_results(counter, output_count_mtx, output_cell_ids, 
                               cell_index, i7_prefix, previous_barcode, gene_exon_id_index)
                # reset
                counter = collections.Counter()
                cell_index += 1
                count_pair_with_junctions(algn_pair, exons, exon_gene_dict, counter)
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
    ligation_barcode, rt_barcode = barcode.split("-")
    output_cell_ids.write(f"{cell_index},{i7_prefix}.{ligation_barcode}.{rt_barcode},{unmatched_rate}\n")

def count_bam_parallel(input_folder, output_folder, exons, gene_exon_id_index, exon_gene_dict, i7_prefix):
    """
    Count exon reads in BAM files for a i7 file prefix (used for parallel processing)
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {i7_prefix} Running...")
    bam_path = os.path.join(input_folder, f"{i7_prefix}.filtered.dedup.bam")
    output_prefix = os.path.join(output_folder, i7_prefix)
    count_bam(bam_path, output_prefix, exons, gene_exon_id_index, exon_gene_dict)
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

    exons, exons_table, exon_gene_dict = get_exons_array(args.gtf)
    out_exon_annotation_path = os.path.join(args.output_folder, "gene_exon_ids.csv.gz" if IS_GZIP else "gene_exon_ids.csv")
    exons_table.to_csv(out_exon_annotation_path, index=False)

    gene_exon_id_index = exons_table.set_index("gene_exon_id")["exon_index"].to_dict()
    del exons_table # free memory

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Exon reference loaded")
    from functools import partial
    from multiprocessing import Pool
    i7_prefixes = list(prefix for prefix in open(args.i7_prefix_file).read().splitlines() if prefix)
    worker_count_bam = partial(count_bam_parallel, args.input_folder, args.output_folder, exons, gene_exon_id_index, exon_gene_dict)
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
    exon_gene_dict = dict() # {gene_exon_id: gene_id}
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
            exon_gene_dict[gene_exon_id] = feature.attr['gene_id']

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

    return exons, exon_count_table, exon_gene_dict

if __name__ == "__main__":
    main()
