#!/usr/bin/env python3

##############################################################
# merging_velocyto.py
# Description: Merge loom files from velocyto into a single loom file.
# This script takes multiple loom files and merges them into one.
#   Parameters:
#     --input: Input loom files (separated by space)
#     --output: Output loom file
##############################################################

import loompy
import argparse

def merge_loom_files(input_files, output_file):
    """
    Merge multiple loom files into a single loom file.
    
    Parameters:
        input_files (list): List of input loom file paths.
        output_file (str): Path to the output loom file.
    """
    loompy.combine(input_files, output_file, key="Accession")

def main():
    parser = argparse.ArgumentParser(description="Merge loom files from velocyto into a single loom file.")
    parser.add_argument("-i", "--input", nargs='+', required=True, help="Input loom files")
    parser.add_argument("-o", "--output", required=True, help="Output loom file")
    
    args = parser.parse_args()
    
    merge_loom_files(args.input, args.output)

    # rename the CellID in loom file
    # e.g. from "Sample1:Cell1" to "Cell1"

    with loompy.connect(args.output, 'r+') as ds:
        ds.ca["CellID"] = [cell.split(":")[1] for cell in ds.ca["CellID"]]

    print(f"Merged {len(args.input)} loom files into {args.output}")

if __name__ == "__main__":
    main()