#!/usr/bin/env python3

##############################################################
# merging_count_gene.py
# Description: Merge gene count files for all i7 demultiplexed files (also merge RT primers)
#   Parameters:
#     --i7_prefix_file <path>: Text file containing i7 demultiplexed file prefixes
#     --input_folder   <path>: Input folder containing count files, _cell_ids, and gene_ids
#     --output_folder  <path>: Output folder to save merged files
#     --rt_barcode_tsv <path>: TSV file containing RT barcodes (Sample	shortdt	randomN)
#     --paired: Use this flag if the pipeline is paired-end
#     --gzip: Use this flag if the count files are gzipped
#   Output:
#     cell_annotation.csv: Cell annotation matrix with metadata (UMI count, gene count, PCR batch, RT ligation barcodes, ligation barcodes, RT barcodes)
#     feature_annotation.csv: Gene annotation matrix
#     expression_matrix.mtx: Sparse expression matrix
##############################################################

# Basicly this script reads the gene count files for all i7 demux files, merges them, 
#   and adds metadata to the cell annotation matrix.
# It also merges RT primers from the same cell
#   (barcodes on the same row are from the same cell)

import pandas as pd
from os import path
from scipy import sparse
import gzip
IS_GZIP = False

def merge_gene_counts(i7_prefixes, input_folder):
    """
    Merge gene count files for i7 demultiplexed files.
    """
    gene_ids_df = pd.read_csv(
        path.join(input_folder, "gene_ids.csv.gz" if IS_GZIP else "gene_ids.csv"),
        index_col=0)

    count_mtx = None
    cell_ids_df = None

    for i7_prefix in i7_prefixes:
        current_count_mtx, current_cell_ids_df = read_matrix(
            path.join(input_folder, i7_prefix), gene_ids_df)
        
        if count_mtx is None:
            count_mtx = current_count_mtx
            cell_ids_df = current_cell_ids_df
        else:
            count_mtx = sparse.hstack([count_mtx, current_count_mtx])
            cell_ids_df = pd.concat([cell_ids_df, current_cell_ids_df])

    return count_mtx, cell_ids_df, gene_ids_df

def merge_exon_intron_counts(gene_ids_df, count_mtx):
    """
    Merge exon and intron counts.
    """
    exon_indices   = gene_ids_df[gene_ids_df["feature"] == "exon"].index
    intron_indices = gene_ids_df[gene_ids_df["feature"] == "intron"].index

    count_mtx_exon   = count_mtx[exon_indices, :]
    count_mtx_intron = count_mtx[intron_indices, :]
    count_mtx_combined = count_mtx_exon + count_mtx_intron

    gene_ids_df_exon = gene_ids_df.loc[exon_indices, ["gene_id", "gene_type", "gene_name"]]

    return gene_ids_df_exon, count_mtx_combined

def read_matrix(file_prefix, gene_ids_df):
    """
    Return a sparse count matrix and cell IDs dataframe for a I7 demultiplexed file.
    """
    count_mtx = pd.read_csv(file_prefix + ".count.gz" if IS_GZIP else ".count")
    cell_ids_df = pd.read_csv(file_prefix + "_cell_ids.csv.gz" if IS_GZIP else "_cell_ids.csv", index_col=0)

    sparse_count_mtx = sparse.coo_matrix(
        (count_mtx["count"], (count_mtx["gene_index"], count_mtx["cell_index"])),
        shape=(gene_ids_df.shape[0], cell_ids_df.shape[0]),
        dtype=int
    )

    return sparse_count_mtx, cell_ids_df

def get_metadata_remove_empty(cell_ids_df, count_mtx):
    """
    Add metadata to the cell annotation matrix (umis, genes, pcr_batch, barcodes)
    And remove cells with 0 reads
    """
    # columns of count_mtx are cells, axis = 0
    cell_ids_df["umi_count"] = count_mtx.sum(axis=0).A[0]
    cell_ids_df["gene_count"] = (count_mtx > 0).sum(axis=0).A[0]
    cell_ids_df[["pcr_batch", "rt_ligation_barcodes"]] = cell_ids_df["cell_id"].str.split(pat=".", n=1, expand=True)
    cell_ids_df["ligation_barcodes"] = cell_ids_df["rt_ligation_barcodes"].str[0:10]
    cell_ids_df["rt_barcodes"] = cell_ids_df["rt_ligation_barcodes"].str[10:20]
    
    count_mtx = count_mtx[:, cell_ids_df["umi_count"] > 0]
    cell_ids_df = cell_ids_df[cell_ids_df["umi_count"] > 0].reset_index(drop=True)

    return cell_ids_df, count_mtx

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

def merge_rt_primers(cell_ids_df, count_mtx, rt_barcode_tsv_path):
    """
    Merge RT primers from one cell
    """
    # barcodes on the same row are from the same cell
    # shortdt randomN (convert randomN to short-dT barcodes for cell merging)
    rt_barcode_dict, randomn_set, shortdt_set, rt_dt_sample_dict = get_rt_barcode_dict(rt_barcode_tsv_path)
    cell_ids_df["rt_barcodes_new"] = cell_ids_df["rt_barcodes"].map(rt_barcode_dict)
    cell_ids_df["cell_id_new"] = cell_ids_df["pcr_batch"] + "." + cell_ids_df["ligation_barcodes"] + cell_ids_df["rt_barcodes_new"]

    # record short-dT and randomN UMI separately
    cell_ids_df["shortdt_umi_count"] = 0.0
    mask = cell_ids_df["rt_barcodes"].isin(shortdt_set)
    cell_ids_df.loc[mask, "shortdt_umi_count"] = cell_ids_df.loc[mask, "umi_count"]
    cell_ids_df["randomn_umi_count"] = 0.0
    mask = cell_ids_df["rt_barcodes"].isin(randomn_set)
    cell_ids_df.loc[mask, "randomn_umi_count"] = cell_ids_df.loc[mask, "umi_count"]

    # handle single and double barcode cells differently
    # Keep cells with only one barcode
    mask_single = cell_ids_df.groupby("cell_id_new")["cell_id_new"].transform("size") == 1
    cell_ids_single = cell_ids_df[mask_single]
    count_mtx_single = count_mtx[:, mask_single]

    # sort double barcode cells by cell_id_new and randomn_umi_count
    count_mtx_double = count_mtx[:, ~mask_single]
    cell_ids_double = cell_ids_df[~mask_single].reset_index(drop=True)
    sorted_indices = cell_ids_double.sort_values(["cell_id_new", "randomn_umi_count"]).index.to_numpy()
    cell_ids_double_sorted = cell_ids_double.loc[sorted_indices].reset_index(drop=True)
    count_mtx_double_sorted = count_mtx_double[:, sorted_indices]
    
    # default sort is ascending, so the first one has 0 randomn_umi_count
    #     So the first one is short-dT, the second one is randomN
    cell_ids_double_shortdt = cell_ids_double_sorted.iloc[::2].reset_index(drop=True)
    count_mtx_double_shortdt = count_mtx_double_sorted[:, ::2]
    cell_ids_double_randomn = cell_ids_double_sorted.iloc[1::2].reset_index(drop=True)
    count_mtx_double_randomn = count_mtx_double_sorted[:, 1::2]

    cell_ids_double_shortdt["randomn_umi_count"] = cell_ids_double_randomn["randomn_umi_count"]
    # cell_ids_double_shortdt["shortdt_umi_count"] += cell_ids_double_randomn["shortdt_umi_count"]
    # cell_ids_double_shortdt["randomn_umi_count"] += cell_ids_double_randomn["randomn_umi_count"]

    # Recalculate 1:unmatched_rate, 2:shortdt_umi_count, 3:gene_count
    total_reads_shortdt = 1 / (1 - cell_ids_double_shortdt["unmatched_rate"]) * cell_ids_double_shortdt["umi_count"]
    total_reads_randomn = 1 / (1 - cell_ids_double_randomn["unmatched_rate"]) * cell_ids_double_randomn["umi_count"]
    cell_ids_double_shortdt["unmatched_rate"] = (cell_ids_double_shortdt["unmatched_rate"] * total_reads_shortdt + cell_ids_double_randomn["unmatched_rate"] * total_reads_randomn) / (total_reads_shortdt + total_reads_randomn)
    cell_ids_double_shortdt["umi_count"] += cell_ids_double_randomn["umi_count"]
    cell_ids_double_shortdt["gene_count"] += cell_ids_double_randomn["gene_count"]

    # add short-dT and randomN UMI counts
    count_mtx_double_shortdt += count_mtx_double_randomn

    # merge single and double barcode cells
    cell_ids_df = pd.concat([cell_ids_single, cell_ids_double_shortdt]).reset_index(drop=True)
    count_mtx = sparse.hstack([count_mtx_single, count_mtx_double_shortdt])

    # edit metadata
    cell_ids_df["gene_count"] = (count_mtx > 0).sum(axis=0).A[0]        # re-calculate gene_count
    cell_ids_df["rt_barcodes_shortdt"] = cell_ids_df["rt_barcodes_new"] # rt_barcodes_new must be short-dT
    cell_ids_df["cell_id"] = cell_ids_df["cell_id_new"]
    cell_ids_df["sample"] = cell_ids_df["rt_barcodes_shortdt"].map(rt_dt_sample_dict)

    # drop unnecessary columns
    cell_ids_df = cell_ids_df.drop(columns=["rt_ligation_barcodes", "rt_barcodes_new", "cell_id_new", "rt_barcodes"])

    return cell_ids_df, count_mtx

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_folder", required=True, help="Input folder containing count files, cell_ids, and gene_ids")
    parser.add_argument("-o", "--output_folder", required=True, help="Output folder to save merged files")
    parser.add_argument("-p", "--i7_prefix_file", required=True, help="Text file containing i7 demultiplexed file prefixes")
    parser.add_argument("-t", "--rt_barcode_tsv", required=True, help="TSV file containing RT barcodes (shortdt	randomN)")
    parser.add_argument("--gzip", action="store_true", help="Use this flag if the count files are gzipped")
    args = parser.parse_args()
    global IS_GZIP
    IS_GZIP = args.gzip

    i7_prefixes = [i7_prefix for i7_prefix in open(args.i7_prefix_file).read().splitlines() if i7_prefix]
    count_mtx, cell_ids_df, gene_ids_df = merge_gene_counts(i7_prefixes, args.input_folder)

    count_mtx = count_mtx.tocsr() # convert to csr format

    # merge exon and intron counts
    gene_ids_df, count_mtx = merge_exon_intron_counts(gene_ids_df, count_mtx)

    # remove empty cells and add metadata (UMI count, gene count, PCR batch, barcodes)
    cell_ids_df, count_mtx = get_metadata_remove_empty(cell_ids_df, count_mtx)
    # merge RT primers
    cell_ids_df, count_mtx = merge_rt_primers(cell_ids_df, count_mtx, args.rt_barcode_tsv)

    from os import makedirs
    makedirs(args.output_folder, exist_ok=True)
    cell_ids_df.to_csv(path.join(args.output_folder, "cell_annotation.csv.gz" if IS_GZIP else "cell_annotation.csv"), index=False)
    gene_ids_df.to_csv(path.join(args.output_folder, "feature_annotation.csv.gz" if IS_GZIP else "feature_annotation.csv"), index=False)
    from scipy.io import mmwrite
    mtx_file = gzip.open(path.join(args.output_folder, "expression_matrix.mtx.gz"), "wb") \
            if IS_GZIP else open(path.join(args.output_folder, "expression_matrix.mtx"), "w")
    mmwrite(mtx_file, count_mtx)
    mtx_file.close()

if __name__ == "__main__":
    main()
