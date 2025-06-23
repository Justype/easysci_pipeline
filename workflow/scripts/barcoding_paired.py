#!/usr/bin/env python3

##############################################################
# barcoding_reads_paired.py
# Description: This script attaches ligation and RT barcodes and UMI to FASTQ files
#   Input:
#      --input_prefix:     Input fastq file prefix  (R1, R2, R3) (R2 is actually I5)
#      --output_prefix:    Output fastq file prefix (R1, R2)
#      --ligation_barcode: Ligation barcode pickle file
#      --rt_barcode:       RT barcode pickle file
#      --rt_barcode_tsv:   RT barcode tsv (Sample	shortdt	randomN)
#      --min_length:       Minimum length of trimmed sequence (default: 20)
#      --barcodes_pattern: Barcodes pattern (ligation barcode in R2, UMI in R1, RT barcode in R1, oligo-dT in R1) as comma-separated values (default: 10,10,8,15)
#      --min_rt_trim_length: Minimum length of RT barcode to be trimmed from R2 (default: 5)
#   Output: R1, R2 fastq files with barcodes attached to headers
#      Name Format: @<ligation_barcode>-<RT_barcode>,UMI,original_read_name
##############################################################

import gzip
import pandas as pd
import pickle
import os
import re

NEXTERA_R1 = "CTGTCTCTTATACACAT"

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Attach Ligation and RT Barcodes and UMI to FASTQ files")
    parser.add_argument("-i", "--input_prefix", type=str, required=True, help="Input fastq file prefix")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output fastq file prefix")
    parser.add_argument("-l", "--ligation_barcode", type=str, required=True, help="Ligation barcode pickle file")
    parser.add_argument("-r", "--rt_barcode", type=str, required=True, help="RT barcode pickle file")
    parser.add_argument("-t", "--rt_barcode_tsv", type=str, required=True, help="RT barcode tsv (Sample	shortdt	randomN)")
    parser.add_argument("-b", "--barcodes_pattern", type=str, default="10,10,8,15", help="Barcodes pattern (ligation barcode in R2, UMI in R1, RT barcode in R1, oligo-dT in R1) as comma-separated values (default: 10,10,8,15)")
    parser.add_argument("--min_rt_trim_length", type=int, default=5, help="Minimum length of RT barcode to be trimmed from R2 (default: 5)")
    parser.add_argument("-m", "--min_length", type=int, default=20, help="Minimum length of trimmed sequence")
    parser.add_argument("--log", help="log file")

    args = parser.parse_args()

    R1_path = args.input_prefix + "_R1.fastq.gz" # Read 1 (UMI + RT barcode + cDNA)
    I5_path = args.input_prefix + "_R2.fastq.gz" # Index 5 (Ligation barcode)
    R2_path = args.input_prefix + "_R3.fastq.gz" # Read 2 (cDNA)
    
    file_check(R1_path)
    file_check(I5_path)
    file_check(R2_path)
    file_check(args.ligation_barcode)
    file_check(args.rt_barcode)
    file_check(args.rt_barcode_tsv)

    # Parse barcodes pattern
    try:
        barcodes_pattern = list(map(int, re.split(r",\s*", args.barcodes_pattern)))
        if len(barcodes_pattern) != 4:
            raise ValueError("Barcodes pattern must contain 4 comma-separated values")
        if any(x <= 0 for x in barcodes_pattern):
            raise ValueError("All values in barcodes pattern must be positive integers")
        # Make sure min_rt_trim_length is less than or equal to the RT barcode length
        if args.min_rt_trim_length > barcodes_pattern[2]:
            raise ValueError("min_rt_trim_length must be less than or equal to the RT barcode length")
    except ValueError as e:
        print(f"Error parsing barcodes pattern: {e}")
        exit(1)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

    ligation_barcode = pickle.load(open(args.ligation_barcode, "rb"))
    rt_barcode = pickle.load(open(args.rt_barcode, "rb"))
    randomN_barcode = pd.read_csv(args.rt_barcode_tsv, sep="\t")["randomN"].tolist()

    with gzip.open(R1_path, "rt", encoding="utf-8") as input_R1, \
         gzip.open(I5_path, "rt", encoding="utf-8") as input_I5, \
         gzip.open(R2_path, "rt", encoding="utf-8") as input_R2, \
         gzip.open(args.output_prefix + "_R1.fastq.gz", "wt", encoding="utf-8") as output_R1, \
         gzip.open(args.output_prefix + "_R2.fastq.gz", "wt", encoding="utf-8") as output_R2:
        
        (n_reads, n_barcode, n_pass, n_rt, n_length) = \
            attach_UMI(input_R1, input_R2, input_I5, output_R1, output_R2, ligation_barcode, rt_barcode, randomN_barcode, args.min_length, barcodes_pattern, args.min_rt_trim_length)
    
    base_name = os.path.basename(args.input_prefix)
    print(f"{base_name}: {n_reads} reads, {n_pass} passed, {n_pass/n_reads:.2%} passed, {n_barcode} barcode failed, {n_rt} RT trimmed, {n_length} length failed")

    if args.log:
        os.makedirs(os.path.dirname(args.log), exist_ok=True)
        # If log file exists, remove it before writing
        if os.path.exists(args.log):
            os.remove(args.log)
        
        with open(args.log, "a") as f:
            f.write(f"{n_reads}\treads_total\n")
            f.write(f"{n_reads - n_barcode}\treads_barcode_passed\n")
            f.write(f"{n_pass}\treads_length_passed\n")

def attach_UMI(input_R1, input_R2, input_I5, output_R1, output_R2, ligation_barcode, rt_barcode, randomN_barcode, min_length, barcodes_pattern=[10, 10, 8, 15], min_rt_trim_length=5):
    """
    Attach Ligation and RT Barcodes and UMI to FASTQ files
    1. Reads without ligation and RT barcodes will be discarded
    2. Ligation and RT barcodes will be corrected to the closest expected barcode
    3. Out Read Name: ligation_barcode + RT_barcode , UMI , original read name
    4. If RT barcode is not randomN_barcode, the sequence will be trimmed more to remove oligo-dT region
    5. Remove RT barcode from R2 (Nextera adapter from R1 will be trimmed later by trim_galore)
    6. If trimmed sequence is less than min_length, the read will be discarded

    min_rt_trim_length: Minimum length of RT barcode to be trimmed from R2 (default: 5)
    """

    n_reads = n_barcode = n_pass = n_nextera = n_rt = n_length = 0

    ligation_length, umi_length, rt_length, polyT_length = barcodes_pattern

    while True:
        read1 = [input_R1.readline().strip() for _ in range(4)]
        if not read1[0]: # End of File
            break
        read2 = [input_R2.readline().strip() for _ in range(4)]
        ligation = [input_I5.readline().strip() for _ in range(4)]
        n_reads += 1

        #1. Make sure the ligation and RT barcodes are correct
        if ligation[1][:ligation_length] not in ligation_barcode or read1[1][umi_length:umi_length+rt_length] not in rt_barcode:
            n_barcode += 1
            continue

        #2. Correct the ligation and RT barcodes
        seq_ligation = ligation_barcode[ligation[1][:ligation_length]]
        seq_rt = rt_barcode[read1[1][umi_length:umi_length+rt_length]]
        seq_umi = read1[1][:umi_length]

        #3. Rename the read name
        #         ligation_barcode RT_barcode, UMI, original read name
        read1[0] = f"@{seq_ligation}-{seq_rt},{seq_umi},{read1[0][1:]}"
        read2[0] = f"@{seq_ligation}-{seq_rt},{seq_umi},{read2[0][1:]}"

        #4. Trim R1 if RT barcode is randomN or oligo-dT
        if seq_rt in randomN_barcode:
            read1[1] = read1[1][umi_length + rt_length:] # remove UMI and RT barcode
            read1[3] = read1[3][umi_length + rt_length:]
        else:
            read1[1] = read1[1][umi_length + rt_length + polyT_length:] # remove UMI, RT barcode and oligo-dT
            read1[3] = read1[3][umi_length + rt_length + polyT_length:]

        # There is no need to trim NEXTERA adapter here, as it will be done by trim_galore later
        # #5. Trim Nextera adapter from R1
        # nextera_index = read1[1].find(NEXTERA_R1)
        # # nextera_index = get_R1_nextera_1mismatch_index(read1[1]) # 1 mismatch allowed (using regex)
        # if nextera_index != -1:
        #     read1[1] = read1[1][:nextera_index]
        #     read1[3] = read1[3][:nextera_index]
        #     n_nextera += 1
        
        #5. Trim RT barcode from R2
        seq_rt_rev = reverse_complement(seq_rt)
        rt_index = read2[1].find(seq_rt_rev)
        if rt_index != -1:
            read2[1] = read2[1][rt_index:]
            read2[3] = read2[3][rt_index:]
            n_rt += 1
        else: # try to trim rt barcode at the end of R2
            for i in range(rt_length, min_rt_trim_length - 1, -1):
                if read2[1].endswith(seq_rt_rev[:i]):
                    read2[1] = read2[1][:-i] # remove the RT barcode from the end
                    read2[3] = read2[3][:-i]
                    n_rt += 1
                    break

        #6. Check if the trimmed sequences are long enough
        if len(read1[1]) < min_length or len(read2[1]) < min_length:
            n_length += 1
            continue

        output_R1.write("\n".join(read1) + "\n")
        output_R2.write("\n".join(read2) + "\n")
        n_pass += 1

    return (n_reads, n_barcode, n_pass, n_rt, n_length)


def file_check(file_path):
    if not os.path.exists(file_path):
        print(f"Error: {file_path} does not exist")
        exit(1)
    
def reverse_complement(seq):
    """
    Reverse complement a DNA sequence
    """
    return seq.translate(str.maketrans("ATCGN", "TAGCN"))[::-1]

def generate_1mismatch_set(seq):
    """
    Generate a set of sequences with 1 mismatch
    """
    seq_set = set()
    for i in range(len(seq)):
        for base in "ATCGN":
            if base != seq[i]:
                seq_set.add(seq[:i] + base + seq[i+1:])
    return seq_set

def get_R2_rt_1mismatch_index(read2, rt_barcode):
    """
    Get the index of the RT barcode with 1 mismatch in Read 2 (-1 if not found)
    """
    match = read2.find(reverse_complement(rt_barcode))
    if match != -1:
        return match
    match = re.search("|".join(generate_1mismatch_set(rt_barcode)), read2)
    if match:
        return match.start()
    return -1

if __name__ == "__main__":
    main()

