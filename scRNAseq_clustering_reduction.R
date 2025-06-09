## --------------------------------------------------------------------------
# Script Name: scRNAseq_Clustering_DimensionalityReduction.R
# Author: Zeng Zhiwei
# Date: 2025-04-11
# Description: This script performs clustering and dimensionality reduction (UMAP, t-SNE)
#               for scRNA-seq data, with optional batch effect correction using
#               Harmony or Seurat's IntegrateData.
# --------------------------------------------------------------------------

# Load required R packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(harmony)
  library(clustree)
  library(future)
  library(optparse)
})

# --------------------------------------------------------------------------
# Set up command line arguments
# --------------------------------------------------------------------------
option_list <- list(
  make_option("--input_dir", type = "character", default = NULL,
              help = "Directory containing the normalized and scaled data (required)."),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Output directory for clustering and reduction results (optional)."),
  make_option("--pc_num", type = "integer", default = 20,
              help = "Number of principal components to use [default= %default]."),
  make_option("--batch_correction", type = "character", default = "harmony",
              help = "Batch effect correction method: 'none', 'harmony', or 'seurat_integrate' [default= %default]."),
  make_option("--batch_variable", type = "character", default = "sample",
              help = "Column name in the Seurat object metadata containing batch information [default= %default]."),
  make_option("--max_resolution", type = "numeric", default = 0.6,
              help = "Maximum resolution for clustering [default= %default]."),
  make_option("--resolution_step", type = "numeric", default = 0.1,
              help = "Step size for clustering resolutions [default= %default]."),
  make_option("--n_cores", type = "integer", default = 4,
              help = "Number of cores to use for parallel processing [default= %default].")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$input_dir)) stop("Error: Input directory (--input_dir) is required.")

# Set up directories and parameters
input_dir <- normalizePath(opt$input_dir, mustWork = TRUE)
output_dir <- opt$output_dir
pc_num <- opt$pc_num
batch_correction <- tolower(opt$batch_correction)
batch_variable <- opt$batch_variable
max_resolution <- opt$max_resolution
resolution_step <- opt$resolution_step
n_cores <- opt$n_cores

# Define output directory
if (is.null(opt$output_dir)) {
  output_dir <- file.path("scRNAseq_results", "clustering_reduction")
} else {
  output_dir <- normalizePath(opt$output_dir, mustWork = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load Seurat object
seurat_file <- file.path(opt$input_dir, "filter_normailze_seurat.rds")
if (!file.exists(seurat_file)) stop("Error: Processed Seurat object file not found: ", seurat_file)
seurat_object_raw <- readRDS(seurat_file)
message("Seurat object loaded successfully.")

# Configure parallel processing
options(future.globals.maxSize = 60 * 1024 * 1024 * 1024)
plan("multicore", workers = n_cores)

# Define clustering resolutions
cluster_resolutions <- seq(0.1, max_resolution, resolution_step)

# Perform batch correction (if applicable) and dimensionality reduction
if (batch_correction == "harmony") {
  message("Applying Harmony batch correction using batch variable: ", batch_variable)
  if (!batch_variable %in% colnames(seurat_object_raw@meta.data)) {
    stop("Error: Batch variable '", batch_variable, "' not found in Seurat object metadata.")
  }
  seurat_object <- RunHarmony(seurat_object_raw, reduction = "pca", theta = 4, dims.use = 1:pc_num,
                              group.by.vars = batch_variable, reduction.save = "harmony")
  reduction_method <- "harmony"
} else if (batch_correction == "seurat_integrate") {
  message("Applying Seurat IntegrateData for batch correction using batch variable: ", batch_variable)
  if (!batch_variable %in% colnames(seurat_object_raw@meta.data)) {
    stop("Error: Batch variable '", batch_variable, "' not found in Seurat object metadata.")
  }
  seurat_list <- SplitObject(seurat_object_raw, split.by = batch_variable)
  integration.anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:pc_num)
  seurat_object <- IntegrateData(anchorset = integration.anchors, dims = 1:pc_num)
  DefaultAssay(seurat_object) <- "integrated"
  seurat_object <- ScaleData(seurat_object, verbose = FALSE)
  seurat_object <- RunPCA(seurat_object, npcs = pc_num, verbose = FALSE)
  reduction_method <- "pca" # PCA on the integrated assay
} else {
  message("Proceeding without batch correction.")
  seurat_object <- seurat_object_raw
  reduction_method <- "pca"
  if (!"pca" %in% names(seurat_object@reductions)) {
    message("PCA not found, running PCA...")
    seurat_object <- RunPCA(seurat_object, npcs = pc_num, verbose = FALSE)
  }
}

# Run dimensionality reduction (if not already done after integration) and clustering
if (!"umap" %in% names(seurat_object@reductions)) {
  seurat_object <- RunUMAP(seurat_object, reduction = reduction_method, dims = 1:pc_num)
}
if (!"tsne" %in% names(seurat_object@reductions)) {
  seurat_object <- RunTSNE(seurat_object, reduction = reduction_method, dims = 1:pc_num)
}
seurat_object <- FindNeighbors(seurat_object, reduction = reduction_method, dims = 1:pc_num)
seurat_object <- FindClusters(seurat_object, resolution = cluster_resolutions)

# Reset parallel processing to sequential
plan("sequential")

# Generate and save plots
p_umap <- DimPlot(seurat_object, reduction = "umap", group.by = batch_variable)
p_tsne <- DimPlot(seurat_object, reduction = "tsne", group.by = batch_variable)
ggsave(file.path(output_dir, paste0(batch_correction, "_umap_plot.pdf")), plot = p_umap)
ggsave(file.path(output_dir, paste0(batch_correction, "_tsne_plot.pdf")), plot = p_tsne)

# Generate and save clustering tree plot
p_clustree <- clustree(seurat_object, prefix = paste0(DefaultAssay(seurat_object), "_snn_res."))
ggsave(file.path(output_dir, paste0(batch_correction, "_clustree_plot.pdf")), plot = p_clustree, height = 14, width = 12)

# Save final Seurat object
saveRDS(seurat_object, file.path(output_dir, paste0(batch_correction, "_clustered_reduced_data.rds")))
message("Clustering and dimensionality reduction completed successfully. Results saved in: ", output_dir)