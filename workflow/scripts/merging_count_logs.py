#!/usr/bin/env python3

import pandas as pd
import os

human_output_path = "human/human_stats.csv"
mouse_output_path = "mouse/mouse_stats.csv"
sample_ids = open("metadata/i7_demux_prefix.txt").read().splitlines()

# reads_total,reads_barcode_passed,reads_length_passed,reads_trimmed,reads_mapped,alignments,alignments_qc_filtered,alignments_dedup
if os.path.exists(os.path.dirname(human_output_path)):
    human_stats = pd.DataFrame()
    for sample in sample_ids:
        if os.path.isfile(f"logs/human/{sample}.log"):
            stats = pd.read_csv(f"logs/human/{sample}.log", sep="\t", index_col=1, header=None)
            stats.columns = [sample]
            human_stats = pd.concat([human_stats, stats], axis=1) # column-wise concatenation
    human_stats = human_stats.fillna(-1).astype(int).transpose()

    # add ratio to reads_total
    human_stats["reads_barcode_passed_ratio"]   = round(human_stats["reads_barcode_passed"] / human_stats["reads_total"], 4)
    human_stats["reads_length_passed_ratio"]    = round(human_stats["reads_length_passed"] / human_stats["reads_total"], 4)
    human_stats["reads_trimmed_ratio"]          = round(human_stats["reads_trimmed"] / human_stats["reads_total"], 4)
    human_stats["reads_mapped_ratio"]           = round(human_stats["reads_mapped"] / human_stats["reads_total"], 4)
    human_stats["alignments_qc_filtered_ratio"] = round(human_stats["alignments_qc_filtered"] / human_stats["reads_total"], 4)
    human_stats["alignments_dedup_ratio"]       = round(human_stats["alignments_dedup"] / human_stats["reads_total"], 4)

    # rearrange columns
    human_stats = human_stats[[
        "reads_total", "reads_barcode_passed", "reads_barcode_passed_ratio",
        "reads_length_passed", "reads_length_passed_ratio",
        "reads_trimmed", "reads_trimmed_ratio",
        "reads_mapped", "reads_mapped_ratio",
        "alignments", "alignments_qc_filtered", "alignments_qc_filtered_ratio",
        "alignments_dedup", "alignments_dedup_ratio"
    ]]
    human_stats.to_csv(human_output_path)

if os.path.exists(os.path.dirname(mouse_output_path)):
    mouse_stats = pd.DataFrame()
    for sample in sample_ids:
        if os.path.isfile(f"logs/mouse/{sample}.log"):
            stats = pd.read_csv(f"logs/mouse/{sample}.log", sep="\t", index_col=1, header=None)
            stats.columns = [sample]
            mouse_stats = pd.concat([mouse_stats, stats], axis=1)
    mouse_stats = mouse_stats.fillna(-1).astype(int).transpose()

    # add ratio to reads_total
    mouse_stats["reads_barcode_passed_ratio"] = round(mouse_stats["reads_barcode_passed"] / mouse_stats["reads_total"], 4)
    mouse_stats["reads_length_passed_ratio"] = round(mouse_stats["reads_length_passed"] / mouse_stats["reads_total"], 4)
    mouse_stats["reads_trimmed_ratio"] = round(mouse_stats["reads_trimmed"] / mouse_stats["reads_total"], 4)
    mouse_stats["reads_mapped_ratio"] = round(mouse_stats["reads_mapped"] / mouse_stats["reads_total"], 4)
    mouse_stats["alignments_qc_filtered_ratio"] = round(mouse_stats["alignments_qc_filtered"] / mouse_stats["reads_total"], 4)
    mouse_stats["alignments_dedup_ratio"] = round(mouse_stats["alignments_dedup"] / mouse_stats["reads_total"], 4)
    
    # rearrange columns
    mouse_stats = mouse_stats[[
        "reads_total", "reads_barcode_passed", "reads_barcode_passed_ratio",
        "reads_length_passed", "reads_length_passed_ratio",
        "reads_trimmed", "reads_trimmed_ratio",
        "reads_mapped", "reads_mapped_ratio",
        "alignments", "alignments_qc_filtered", "alignments_qc_filtered_ratio",
        "alignments_dedup", "alignments_dedup_ratio"
    ]]
    mouse_stats.to_csv(mouse_output_path)
