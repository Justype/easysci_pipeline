#!/usr/bin/env python3

# Read sam from STDIN and write sam to STDOUT
# Usage samtools view -h file.bam | python3 duplicate_removal_paired.py | samtools view -bS - > file_filtered.bam
#   If you want to log the number of reads, use --log option

# The SAM are expected to be sorted by read name
# Name Format: @<ligation_barcode><RT_barcode>,UMI,original_read_name
#      So the reads with the same barcode and UMI are expected to be adjacent

# Only alignment with the same barcode, UMI, and start position of the read pairs are considered duplicates

import argparse
import sys
import os
import pysam

def main():
    parser = argparse.ArgumentParser(description="Remove duplicates based on the cell identity, UMI sequence and mapping locations")
    parser.add_argument("--log", help="log file")
    args = parser.parse_args()

    in_bam = pysam.AlignmentFile("-", "rb") # read bam from stdin
    out_bam = pysam.AlignmentFile("-", "wb", template=in_bam) # write bam to stdout

    pre_locations = set() # locations of the previous read pair
    pre_barcode_UMI = ""  # the barcode and UMI of the previous read pair
    n_algn = n_pass = 0

    while True:
        try:
            algn1 = next(in_bam)
            algn2 = next(in_bam)
        except StopIteration:
            break

        n_algn += 1 # 1 pair of reads

        if algn1.query_name != algn2.query_name:
            sys.stderr.write("WARNING: not proper read pairing during duplicate removal\n")
            continue

        current_barcode_UMI = ",".join(algn1.query_name.partition(",")[:2])
        current_location = f"{algn1.reference_name}:{algn1.reference_start}-{algn2.reference_name}:{algn2.reference_start}"

        if current_barcode_UMI == pre_barcode_UMI: # same barcode and UMI (may be duplicates)
            if current_location in pre_locations:
                continue
            else:
                pre_locations.add(current_location)
        else: # different barcode and UMI
            pre_barcode_UMI = current_barcode_UMI   # reset the barcode and UMI
            pre_locations = set([current_location])
        
        n_pass += 1
        out_bam.write(algn1)
        out_bam.write(algn2)

    in_bam.close()
    out_bam.close()

    sys.stderr.write(f"\r{n_algn} reads, {n_pass} passed, {n_pass/n_algn:.2%} passed\n")
    if args.log:
        os.makedirs(os.path.dirname(args.log), exist_ok=True)
        with open(args.log, "a") as f: # append stat to the log file
            f.write(f"{n_algn}\talignments_qc_filtered\n")
            f.write(f"{n_pass}\talignments_dedup\n")

if __name__ == "__main__":
    main()
