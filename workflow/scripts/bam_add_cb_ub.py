# Add cell barcodes (CB) and UMI (UB) tags to BAM files
# accept BAM from stdin and write to stdout
# PCR barcodes are provided by passing -p --i7_prefix

import argparse
import sys
import os
import pysam
import pandas as pd

def get_rt_barcode_dict(rt_barcode_tsv_path):
    """
    Return a dictionary mapping randomN to short-dT barcodes
    """
    mapping_dict = {}
    randomn_set = set()
    shortdt_set = set()

    rt_barcode_df = pd.read_csv(rt_barcode_tsv_path, sep="\t")

    for _, row in rt_barcode_df.iterrows():
        randomn_set.add(row["randomN"])
        shortdt_set.add(row["shortdt"])
        mapping_dict[row["randomN"]] = row["shortdt"] # randomN -> short-dT
        mapping_dict[row["shortdt"]] = row["shortdt"] # keep short-dT -> short-dT

    return mapping_dict, randomn_set, shortdt_set, rt_barcode_df.set_index("shortdt")["Sample"].to_dict()

def main():
    parser = argparse.ArgumentParser(description="Add cell barcodes (CB) and UMI (UB) tags to BAM files")
    parser.add_argument("-p", "--i7_prefix", help="PCR barcode prefix, e.g. 'ATCACG'", required=True)
    parser.add_argument("-t", "--rt_barcode_tsv", help="Path to the RT barcode TSV file", required=True)
    args = parser.parse_args()

    if not os.path.exists(args.rt_barcode_tsv):
        print(f"RT barcode TSV file {args.rt_barcode_tsv} does not exist. Please check the path.")
        sys.exit(1)
    
    # Get the mapping dictionary and sets from the RT barcode TSV file
    randomn_dt_dict, randomn_set, shortdt_set, sample_dict = get_rt_barcode_dict(args.rt_barcode_tsv)

    in_bam = pysam.AlignmentFile("-", "rb")  # read bam from stdin
    out_bam = pysam.AlignmentFile("-", "wb", template=in_bam)  # write bam to stdout

    for algn in in_bam:
        algn_barcode, algn_UMI = algn.query_name.split(",")[:2]
        # Algn Name Format: <ligation_barcode>-<RT_barcode>,UMI,original_read_name
        algn_ligation, algn_rt_barcode = algn_barcode.split("-")
        algn_rt_barcode = randomn_dt_dict.get(algn_rt_barcode, None)  # convert randomN to short-dT
        if algn_rt_barcode is None:
            print(f"Warning: RT barcode {algn_rt_barcode} not found in the mapping. Skipping alignment.")
            continue
        algn_barcode = f"{args.i7_prefix}-{algn_ligation}-{algn_rt_barcode}"
        algn.set_tag("CB", algn_barcode)
        algn.set_tag("UB", algn_UMI)

        out_bam.write(algn)

    in_bam.close()
    out_bam.close()

if __name__ == "__main__":
    main()
