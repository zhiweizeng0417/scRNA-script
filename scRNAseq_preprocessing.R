# --------------------------------------------------------------------------
# Script Name: scRNAseq_Preprocessing.R
# Author: Zeng Zhiwei
# Date: 2025-04-10
# Description: This script performs preprocessing for single-cell RNA-seq data.
# --------------------------------------------------------------------------

# Load required libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(optparse)

# --------------------------------------------------------------------------
# Define command-line arguments
# --------------------------------------------------------------------------
option_list <- list(
  make_option(c("-n", "--sample_names"), type = "character", default = NULL,
              help = "Comma-separated list of sample names (required)."),
  make_option(c("--matrix_dir"), type = "character", default = NULL,
              help = "Directory containing matrix data for each sample (required)."),
  make_option(c("--doublet_dir"), type = "character", default = NULL,
              help = "Directory containing doublet files for each sample (required)."),
  make_option(c("-o", "--output_dir"), type = "character", default = "scRNAseq_Preprocess",
              help = "Output directory for results."),
  make_option(c("--mt_pattern"), type = "character", default = NULL,
              help = "Pattern for mitochondrial gene names (optional)."),
  make_option(c("--pt_pattern"), type = "character", default = NULL,
              help = "Pattern for chloroplast gene names (optional).")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
# --------------------------------------------------------------------------
validate_parameters <- function(opt) {
  # Validate that required arguments are provided
  required_args <- c("sample_names", "matrix_dir", "doublet_dir")
  missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]
  
  if (length(missing_args) > 0) {
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
  }
  
  # Validate that file directory are exist
  validate_directory <- function(path, name) {
    if (!dir.exists(path)) {
      stop(sprintf("Directory does not exist for --%s: %s", name, path))
    }
  }
  validate_directory(opt$matrix_dir, "matrix_dir")
  validate_directory(opt$doublet_dir, "doublet_dir")
}

validate_parameters(opt)

# Define paths and parameters
sample_names <- unlist(strsplit(opt$sample_names, ","))
matrix_dir <- opt$matrix_dir
doublet_dir <- opt$doublet_dir
output_dir <- opt$output_dir
mt_pattern <- opt$mt_pattern
pt_pattern <- opt$pt_pattern

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# Process each sample
# --------------------------------------------------------------------------
seurat_list <- list()
metadata <- c()

for (sample in sample_names) {
  sample_matrix_path <- file.path(matrix_dir, sample)
  sample_doublet_path <- file.path(doublet_dir, paste0(sample, "_doublets.txt"))
  
  if (!dir.exists(sample_matrix_path)) stop(paste("Matrix directory missing:", sample_matrix_path))
  if (!file.exists(sample_doublet_path)) stop(paste("Doublet file missing:", sample_doublet_path))
  
  seurat_list[[sample]] <- Read10X(data.dir = sample_matrix_path)
  total_cells_before <- ncol(seurat_list[[sample]])
  
  doublet_info <- read.table(sample_doublet_path, header = TRUE, stringsAsFactors = FALSE)
  doublet_barcodes <- doublet_info$barcode[doublet_info$scrublet_DropletType == "doublet"]
  
  seurat_list[[sample]] <- seurat_list[[sample]][, !colnames(seurat_list[[sample]]) %in% doublet_barcodes]
  total_cells_after <- ncol(seurat_list[[sample]])
  
  metadata <- c(metadata, rep(sample, total_cells_after))
  
  write.table(data.frame(Sample = sample, Before = total_cells_before, After = total_cells_after),
              file = file.path(output_dir, paste0(sample, "_cell_counts.txt")),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  colnames(seurat_list[[sample]]) <- paste0(sample, "_", colnames(seurat_list[[sample]]))
}

# --------------------------------------------------------------------------
# Merge Seurat objects
# --------------------------------------------------------------------------
if (length(seurat_list) > 1) {
  combined_matrix <- do.call(cbind, seurat_list)
  sample_mapping <- unlist(lapply(names(seurat_list), function(s) rep(s, ncol(seurat_list[[s]]))))
} else {
  combined_matrix <- seurat_list[[1]]
  sample_mapping <- rep(names(seurat_list)[1], ncol(seurat_list[[1]]))
}

metadata_df <- data.frame(Sample = sample_mapping, row.names = colnames(combined_matrix))
seurat_object <- CreateSeuratObject(counts = combined_matrix, project = "scRNAseq", meta.data = metadata_df, 
                                    min.cells = 10, min.features = 200)

# --------------------------------------------------------------------------
# QC Metrics & Visualization
# --------------------------------------------------------------------------
seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)

if (!is.null(mt_pattern)) seurat_object[['percent.mt']] <- PercentageFeatureSet(seurat_object, pattern = mt_pattern)
if (!is.null(pt_pattern)) seurat_object[['percent.pt']] <- PercentageFeatureSet(seurat_object, pattern = pt_pattern)

features_to_plot <- c('nFeature_RNA', 'nCount_RNA', 'log10GenesPerUMI')
if (!is.null(mt_pattern)) features_to_plot <- c(features_to_plot, 'percent.mt')
if (!is.null(pt_pattern)) features_to_plot <- c(features_to_plot, 'percent.pt')

p <- VlnPlot(seurat_object, features_to_plot, group.by = 'Sample')
ggsave(filename = file.path(output_dir, 'QC_Violin_Plots.pdf'))

# --------------------------------------------------------------------------
# Save the Seurat object
# --------------------------------------------------------------------------
saveRDS(seurat_object, file = file.path(output_dir, 'preprocess_seurat.rds'))
print("Preprocessing completed successfully.")
