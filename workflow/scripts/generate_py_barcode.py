#!/bin/env python3

# required packages: pandas openpyxl
import argparse
from os import path
import pickle
import pandas as pd 

# INPUT: barcode.xlsx (Sheet: Ligation, randomN, dT)

# I5 will be reverse complement of ligation barcode
# R1 will begin with 8bp UMI, followed by 10bp RT barcode
# RT barcode is also a Well barcode (Each Well can be human or mouse sample)

# What will this script do?
# 1. Generate 1-mismatch barcode dictionary for ligation barcodes
# 2. Generate 1-mismatch barcode dictionary for RT barcodes
# 3. Generate human and mouse RT barcode tables (shortdt	randomN)

def main():
    parser = argparse.ArgumentParser(description="Generate 1-mismatch barcode dictionary for ligation and RT barcodes")
    parser.add_argument("barcode_xlsx", help="Path to barcode xlsx file")
    parser.add_argument("output_dir", help="Path to output directory")
    parser.add_argument("--forbid-na", action="store_true", help="Forbid NA values in Sample column")

    args = parser.parse_args()
    df_barcodes = pd.read_excel(args.barcode_xlsx, sheet_name=None)
    # base_dir = path.dirname(args.barcode_xlsx) # Use the directory of the barcode xlsx file as base directory
    base_dir = args.output_dir # Use the output directory as base directory
    output_ligation       = path.join(base_dir, "ligation_barcodes.pkl")
    output_RT_human       = path.join(base_dir, "RT_barcodes_human.pkl")
    output_RT_mouse       = path.join(base_dir, "RT_barcodes_mouse.pkl")
    output_RT_human_table = path.join(base_dir, "RT_barcodes_human.tsv") # shortdt	randomN
    output_RT_mouse_table = path.join(base_dir, "RT_barcodes_mouse.tsv") # shortdt	randomN

    #region Ligation
    ligation_list = df_barcodes["Ligation"]["Barcode"].tolist()
    ligation_revcomp_list = [reverse_complement(barcode) for barcode in ligation_list] # Do not forget to reverse complement
    ligation_dict = generate_barcode_dict(ligation_revcomp_list)
    pickle.dump(ligation_dict, open(output_ligation, "wb"))
    print("============== Ligation Barcodes ==============")
    print(f"Ligation barcodes:   {len(ligation_revcomp_list)}")
    print(f"Ligation 1-mismatch: {len(ligation_dict)}")
    print(f"Ratio:               {len(ligation_dict) / len(ligation_revcomp_list):.1f}")
    print(f"Expected Ratio:      {1 + len(ligation_revcomp_list[0]) * 4.0} = 1 + len ({len(ligation_revcomp_list[0])}) * 4")
    if len(ligation_dict) / len(ligation_revcomp_list) != 1 + len(ligation_revcomp_list[0]) * 4.0:
        print("ERROR: Ligation 1-mismatch ratio is not as expected")
        exit(1)
    #endregion

    #region RT
    # .copy() to convert slice to dataframe
    df_randomN = df_barcodes["randomN"][["Well", "Sample", "Barcode"]].copy()
    df_randomN.rename(columns={"Barcode": "randomN"}, inplace=True) # rename Barcode to randomN

    df_dT = df_barcodes["dT"][["Well", "Barcode"]].copy()
    df_dT.rename(columns={"Barcode": "shortdt"}, inplace=True)  # rename Barcode to dT

    # Join the two dataframes
    df_rt_merge = pd.merge(df_randomN, df_dT, on="Well")
    del df_randomN, df_dT
    df_rt_merge.drop("Well", axis=1, inplace=True) # drop Well column

    if not args.forbid_na:
        df_rt_merge.dropna(inplace=True) # drop rows with NaN values
    
    df_rt_merge["Human"] = df_rt_merge["Sample"].apply(lambda x: x.split("_")[1])
    df_rt_merge["Mouse"] = df_rt_merge["Sample"].apply(lambda x: x.split("_")[2])
    df_rt_merge["Sample"] = df_rt_merge["Sample"].apply(lambda x: "_".join(x.split("_")[3:]))

    print("========= Human and Mouse RT Barcodes =========")

    rt_human_df = df_rt_merge[df_rt_merge["Human"] == "1"][["Sample", "shortdt", "randomN"]]
    rt_human_list = rt_human_df["shortdt"].tolist() + rt_human_df["randomN"].tolist()
    if rt_human_list:
        rt_human_dict = generate_barcode_dict(rt_human_list)
        print(f"Human RT barcodes: {len(rt_human_list)} = {len(rt_human_df)} * 2 (shortdt, randomN)")
        print(f"Human 1-mismatch:  {len(rt_human_dict)}")
        print(f"Ratio:             {len(rt_human_dict) / len(rt_human_list):.1f}")
        print(f"Expected Ratio:    {1 + len(rt_human_list[0]) * 4.0} = 1 + len ({len(rt_human_list[0])}) * 4")
        if len(rt_human_dict) / len(rt_human_list) != 1 + len(rt_human_list[0]) * 4.0:
            print("ERROR: Human 1-mismatch ratio is not as expected")
            exit(1)

        rt_human_df.to_csv(output_RT_human_table, sep="\t", index=False)
        pickle.dump(rt_human_dict, open(output_RT_human, "wb"))
        print("===============================================")

    rt_mouse_df = df_rt_merge[df_rt_merge["Mouse"] == "1"][["Sample", "shortdt", "randomN"]]
    rt_mouse_list = rt_mouse_df["shortdt"].tolist() + rt_mouse_df["randomN"].tolist()
    if rt_mouse_list:
        rt_mouse_dict = generate_barcode_dict(rt_mouse_list)
        print(f"Mouse RT barcodes: {len(rt_mouse_list)} = {len(rt_mouse_df)} * 2 (shortdt, randomN)")
        print(f"Mouse 1-mismatch:  {len(rt_mouse_dict)}")
        print(f"Ratio:             {len(rt_mouse_dict) / len(rt_mouse_list):.1f}")
        print(f"Expected Ratio:    {1 + len(rt_mouse_list[0]) * 4.0} = 1 + len ({len(rt_mouse_list[0])}) * 4")
        if len(rt_mouse_dict) / len(rt_mouse_list) != 1 + len(rt_mouse_list[0]) * 4.0:
            print("ERROR: Mouse 1-mismatch ratio is not as expected")
            exit(1)
        
        rt_mouse_df.to_csv(output_RT_mouse_table, sep="\t", index=False)
        pickle.dump(rt_mouse_dict, open(output_RT_mouse, "wb"))
    #endregion

#region Functions
def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of a DNA sequence
    """
    return seq.translate(str.maketrans("ATCGN", "TAGCN"))[::-1]

def generate_1_mismatch_set(seq: str) -> set:
    """
    Generate a set of 1-mismatch sequences for a given sequence
    """
    mismatch = set([seq])
    for i in range(len(seq)):
        for base in "ATCGN":
            if base != seq[i]:
                mismatch.add(seq[:i] + base + seq[i+1:])
    return mismatch

def generate_barcode_dict(barcode_list: list|set) -> dict:
    """
    Generate { 1-mismatch : original barcode } dictionary
    """
    return {
        barcode_1mismatch : barcode
        for barcode in barcode_list
            for barcode_1mismatch in generate_1_mismatch_set(barcode)
    }
#endregion

if __name__ == "__main__":
    main()

