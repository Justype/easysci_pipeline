#!/usr/bin/env Rscript

library(tidyverse, quietly = TRUE)
library(Seurat, quietly = TRUE)

# Read from the command line
# arg1 input folder
# arg2 output rds file
# arg3 gzip or not

args <- commandArgs(trailingOnly = TRUE)

input_folder <- args[1]
output_file <- args[2]
gzip <- ifelse(length(args) > 2, args[3] == "gzip", FALSE)

cell_path <- ifelse(gzip, "cell_annotation.csv.gz", "cell_annotation.csv")
cell_path <- file.path(input_folder, cell_path)
feature_path <- ifelse(gzip, "feature_annotation.csv.gz", "feature_annotation.csv")
feature_path <- file.path(input_folder, feature_path)
mtx_path <- ifelse(gzip, "expression_matrix.mtx.gz", "expression_matrix.mtx")
mtx_path <- file.path(input_folder, mtx_path)

cells <- read_csv(cell_path, show_col_types = FALSE)

mtx <- ReadMtx(
  mtx = mtx_path, features = feature_path, cells = cell_path,
  feature.sep = ",", feature.column = 3, skip.feature = 1, # Column 3 is gene name
  cell.sep = ",", cell.column = 1, skip.cell = 1
)

so <- CreateSeuratObject(counts = mtx)

# Add metadata from the cell annotation file
so@meta.data %>%
  rownames_to_column("cell_id") %>%
  left_join(cells, by = join_by("cell_id")) %>%
  column_to_rownames("cell_id") -> so@meta.data

saveRDS(so, output_file)
