#!/usr/bin/env Rscript

library(tidyverse, quietly = TRUE)
library(Seurat, quietly = TRUE)

# Read from the command line
# arg1 input folder
# arg2 output rds file
# arg3 gzip or not (by default TRUE)

args <- commandArgs(trailingOnly = TRUE)

input_folder <- args[1]
# Remove trailing slash if it exists
input_folder <- str_remove(input_folder, "/$")
output_file <- ifelse(length(args) > 1, args[2], file.path(dirname(input_folder), paste0(basename(input_folder), ".qs")))
message("Output file: ", output_file)
gzip <- ifelse(length(args) > 2, args[3] != "gzip", TRUE)

read_count_matrix <- function(folder, is_gzipped = TRUE) {
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

feature_annotation <- read_csv(
  file.path(input_folder, ifelse(gzip, "feature_annotation.csv.gz", "feature_annotation.csv")),
  show_col_types = FALSE
)

is_human <- str_detect(feature_annotation$gene_id[1], "^ENSG")
is_mouse <- str_detect(feature_annotation$gene_id[1], "^ENSMUSG")

if (!is_human && !is_mouse) {
  stop("The gene IDs in the feature annotation file must be either human (ENSG) or mouse (ENSMUSG).")
} else if (is_human) {
  message("Detected human gene IDs (ENSG).")
} else if (is_mouse) {
  message("Detected mouse gene IDs (ENSMUSG).")
}

mt_gene_ids <- feature_annotation %>%
  filter(str_detect(gene_name, ifelse(is_human, "^MT-", "^mt-"))) %>%
  pull(gene_id)

message("[", date(), "] Reading count matrix from folder: ", input_folder)

so <- read_count_matrix(input_folder, gzip)

message("[", date(), "] Read successfully. Ncells: ", ncol(so), ", Nfeatures: ", nrow(so))

so[["percent.mt"]] <- PercentageFeatureSet(so, features = mt_gene_ids)

# extension qs or rds
if (str_detect(output_file, regex("\\.qs$", ignore_case = TRUE))) {
  qs::qsave(so, output_file)
} else if (str_detect(output_file, regex("\\.rds$", ignore_case = TRUE))) {
  saveRDS(so, output_file)
} else {
  stop("Output file must have .qs or .rds extension")
}

message("[", date(), "] Saved Seurat object to: ", output_file)
