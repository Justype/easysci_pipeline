#!/bin/bash
# NOTE: This is just an example, do not use it as-is

# Run Snakemake with conda environment
eval "$(conda shell.bash hook)"
conda activate sci_rna
snakemake \
    --cores 96 \
    --rerun-triggers mtime

# Must have pandas installed
python workflow/utils/merging_count_logs.py

# Packages: Seurat tidyverse qs
conda activate r_seurat
# Rscript workflow/utils/read_gene_exon.R output/final/human
#                                        input_prefix        output qs file    min_feature_rna
Rscript workflow/utils/read_gene_exon.R output/final/human output/final/human.qs 100

exit

# Run with Snakemake modules 
module purge
# Specify the version or leave it blank for the default version
module load snakemake

snakemake \
    --use-conda \
    --cores 96 \
    --rerun-triggers mtime
