#!/bin/bash

# Run Snakemake with conda environment
eval "$(conda shell.bash hook)"
conda activate sci_rna
snakemake \
    --cores 96 \
    --rerun-triggers mtime

# Must have pandas installed
python workflow/utils/merging_count_logs.py

exit

# Run with Snakemake modules 
module purge
# Specify the version or leave it blank for the default version
module load snakemake

snakemake \
    --use-conda \
    --cores 96 \
    --rerun-triggers mtime
