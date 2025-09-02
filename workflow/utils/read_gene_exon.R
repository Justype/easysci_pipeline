# Functions ----
.generate_easysci_paths_vector <- function(folder, is_gzipped = TRUE) {
  folder <- sub("/$", "", folder) # remove folder slash if exists
  cell_path <- ifelse(is_gzipped, "cell_annotation.csv.gz", "cell_annotation.csv")
  cell_path <- file.path(folder, cell_path)
  feature_path <- ifelse(is_gzipped, "feature_annotation.csv.gz", "feature_annotation.csv")
  feature_path <- file.path(folder, feature_path)
  mtx_path <- ifelse(is_gzipped, "expression_matrix.mtx.gz", "expression_matrix.mtx")
  mtx_path <- file.path(folder, mtx_path)

  return(c(
    cell_path = cell_path,
    feature_path = feature_path,
    mtx_path = mtx_path
  ))
}

read_gene_exon <- function(prefix = "final/human", gene_folder_suffix = "_gene", exon_folder_suffix = "_exon", is_gzipped = TRUE, symbol = FALSE, min_feature_rna = 0, verbose = TRUE) {
  require(Seurat); require(readr); require(tidyr); require(dplyr); require(stringr)

  t_start <- Sys.time()
  prefix <- sub("/$", "", prefix) # remove folder slash if exists
  gene_paths <- .generate_easysci_paths_vector(paste0(prefix, gene_folder_suffix), is_gzipped)
  exon_paths <- .generate_easysci_paths_vector(paste0(prefix, exon_folder_suffix), is_gzipped)

  gene_features <- read_csv(gene_paths["feature_path"], show_col_types = FALSE)
  is_human <- str_detect(gene_features$gene_id[1], "^ENSG")
  is_mouse <- str_detect(gene_features$gene_id[1], "^ENSMUSG")
  if (!is_human && !is_mouse) {
    stop("The gene IDs in the feature annotation file must be either human (ENSG) or mouse (ENSMUSG).")
  } else if (is_human) {
    if (verbose) message("Detected human gene IDs (ENSG).")
  } else if (is_mouse) {
    if (verbose) message("Detected mouse gene IDs (ENSMUSG).")
  }

  cells <- read_csv(gene_paths["cell_path"], show_col_types = FALSE)

  # get gene and exon count matrix
  if (verbose) message("Loading gene matrix...")
  t0 <- Sys.time()
  gene_mtx <- ReadMtx(
    mtx = gene_paths["mtx_path"], features = gene_paths["feature_path"], cells = gene_paths["cell_path"],
    feature.sep = ",", feature.column = ifelse(symbol, 3, 1), skip.feature = 1, # 1 ENSEMBL ID; 3 GENE SYMBOL
    cell.sep = ",", cell.column = 1, skip.cell = 1
  )
  if (verbose) message("Gene matrix loaded in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " minutes.")

  if (verbose) message("Loading exon matrix...")
  t1 <- Sys.time()
  exon_mtx <- ReadMtx(
    mtx = exon_paths["mtx_path"], features = exon_paths["feature_path"], cells = exon_paths["cell_path"],
    feature.sep = ",", feature.column = ifelse(symbol, 3, 1), skip.feature = 1, # 1 <ENSEMBL GENE ID>-<EXON ID>; 3 <GENE SYMBOL>-<ENSEMBL EXON ID>
    cell.sep = ",", cell.column = 1, skip.cell = 1
  )
  if (verbose) message("Exon matrix loaded in ", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 2), " minutes.")

  common_cells <- intersect(colnames(gene_mtx), colnames(exon_mtx))
  gene_mtx <- gene_mtx[, common_cells]
  exon_mtx <- exon_mtx[, common_cells]

  if (min_feature_rna > 0) {
    if (verbose) message("Filtering cells with fewer than ", min_feature_rna, " genes...")
    keep_cells <- colSums(gene_mtx > 0) >= min_feature_rna
    gene_mtx <- gene_mtx[, keep_cells]
    exon_mtx <- exon_mtx[, keep_cells] # keep same cells in exon assay
    if (verbose) message(sum(keep_cells), " / ", length(keep_cells), " cells remaining after filtering.")
  }

  if (verbose) message("Creating Seurat object with cell info and percent.mt...")
  # Add exon counts as a new assay
  seurat_obj <- CreateSeuratObject(counts = gene_mtx, assay = "RNA")
  seurat_obj[["EXON"]] <- CreateAssayObject(counts = exon_mtx)

  # Add metadata
  seurat_obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    left_join(cells, by = join_by("cell_id")) %>%
    column_to_rownames("cell_id") -> seurat_obj@meta.data

  mt_pattern <- ifelse(is_human, "^MT-", "^mt-")
  if (symbol) {
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  } else {
    mt_gene_ids <- gene_features %>%
      filter(str_detect(gene_name, mt_pattern)) %>%
      pull(gene_id)

    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mt_gene_ids)
  }
  if (verbose) {
    total_time <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 2)
    message("Finished creating Seurat object with RNA and EXON assays in ", total_time, " minutes.")
    message("RNA assay: ", nrow(seurat_obj[["RNA"]]), " genes x ", ncol(seurat_obj[["RNA"]]), " cells")
    message("EXON assay: ", nrow(seurat_obj[["EXON"]]), " exons x ", ncol(seurat_obj[["EXON"]]), " cells")
  }

  return(seurat_obj)
}

# Read arguments from command line ----
library(tidyverse)
library(Seurat)

# Read from the command line
# arg1 input prefix
# arg2 output rds/qs file (optional)
# arg3 minimum number of features (optional)

args <- commandArgs(trailingOnly = TRUE)
default_min_feature_rna <- 10

if (length(args) < 1) {
  stop("At least one argument (input prefix) must be provided.")
}
input_prefix <- args[1]
# Remove trailing slash if it exists
input_prefix <- str_remove(input_prefix, "/$")
output_file <- ifelse(length(args) > 1, args[2], file.path(dirname(input_prefix), paste0(basename(input_prefix), ".qs")))
min_feature_rna <- ifelse(length(args) > 2, as.numeric(args[3]), default_min_feature_rna)
message("Output file: ", output_file)

so <- read_gene_exon(prefix = input_prefix, min_feature_rna = min_feature_rna, verbose = TRUE)

t_save <- Sys.time()
message("Saving Seurat object to: ", output_file)
if (str_detect(output_file, regex("\\.qs$", ignore_case = TRUE))) {
  qs::qsave(so, output_file)
} else if (str_detect(output_file, regex("\\.rds$", ignore_case = TRUE))) {
  saveRDS(so, output_file)
} else {
  stop("Output file must have .qs or .rds extension")
}

message("Seurat object saved to: ", output_file, " in ", round(as.numeric(difftime(Sys.time(), t_save, units = "mins")), 2), " minutes.")
