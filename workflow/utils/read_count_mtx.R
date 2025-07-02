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

read_count_matrix <- function(folder, is_gzipped = T) {
  require(Seurat); require(readr); require(tidyr); require(dplyr)
  cell_path <- ifelse(is_gzipped, "cell_annotation.csv.gz", "cell_annotation.csv")
  cell_path <- file.path(folder, cell_path)
  feature_path <- ifelse(is_gzipped, "feature_annotation.csv.gz", "feature_annotation.csv")
  feature_path <- file.path(folder, feature_path)
  mtx_path <- ifelse(is_gzipped, "expression_matrix.mtx.gz", "expression_matrix.mtx")
  mtx_path <- file.path(folder, mtx_path)
  
  cells <- read_csv(cell_path, show_col_types = FALSE)
  
  mtx <- ReadMtx(
    mtx = mtx_path, features = feature_path, cells = cell_path,
    feature.sep = ",", feature.column = 1, skip.feature = 1, # 1 ENSEMBL ID; 3 GENE SYMBOL
    cell.sep = ",", cell.column = 1, skip.cell = 1
  )
  
  so <- CreateSeuratObject(counts = mtx)
  
  # Add metadata from the cell annotation file
  so@meta.data %>%
    rownames_to_column("cell_id") %>%
    left_join(cells, by = join_by("cell_id")) %>%
    column_to_rownames("cell_id") -> so@meta.data
  
  return(so)
}

saveRDS(read_count_matrix(input_folder, gzip), output_file)
