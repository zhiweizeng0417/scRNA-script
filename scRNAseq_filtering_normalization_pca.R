#--------------------------------------------------------------------------
# Script Name: scRNAseq_filtering_normalizztion_pca.R
# Author: ZengZhiwei
# Date: 2025-04-10
# Description: This script performs cell filtering, normalization, and PCA for scRNA-seq data.
# --------------------------------------------------------------------------

# Load required R packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(optparse)
})

# --------------------------------------------------------------------------
# Set up command line arguments
# --------------------------------------------------------------------------
option_list <- list(
  make_option("--input_dir", type = "character", default = NULL,
              help = "Directory containing the Seurat object RDS file (required)."),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Output directory [default: 'scRNAseq_results/filter_normailze']"),
  make_option("--lower_quantile", type = "numeric", default = 0.05,
              help = "Lower percentile for nFeature_RNA and nCount_RNA filtering [default= %default]."),
  make_option("--upper_quantile", type = "numeric", default = 0.95,
              help = "Upper percentile for nFeature_RNA filtering [default= %default]."),
  make_option("--max_mt_percent", type = "numeric", default = 20,
              help = "Maximum allowed percentage of mitochondrial genes [default= %default]."),
  make_option("--max_pt_percent", type = "numeric", default = 20,
              help = "Maximum allowed percentage of chloroplast genes [default= %default]."),
  make_option("--normalization_method", type = "character", default = "SCT",
              help = "Normalization method: 'SCT' or 'LogNormalize' [default= %default]."),
  make_option("--nfeatures", type = "integer", default = 2000,
              help = "Number of features to select for variable features [default= %default]."),
  make_option("--log10GenesPerUMI_threshold", type = "numeric" ,default = 0.85,
              help = "The threshold of log10GenesPerUMI(Library Complexity) [default= %default].")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --------------------------------------------------------------------------
# Validate required arguments
# --------------------------------------------------------------------------
validate_parameters <- function(opt) {
  if (is.null(opt$input_dir)) stop("Input directory (--input_dir) is required")
  if (!dir.exists(opt$input_dir)) stop("Input directory not found: ", opt$input_dir)
  if (opt$lower_quantile < 0 | opt$upper_quantile > 1) stop("Quantiles must be between 0 and 1")
  if (!opt$normalization_method %in% c("SCT", "LogNormalize")) {
    stop("Normalization method must be SCT or LogNormalize")
  }
}
validate_parameters(opt)

# Set up output directory
output_dir <- ifelse(is.null(opt$output_dir),
                     file.path("scRNAseq_results", "filter_normalize"),
                     opt$output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Enable parallel processing
plan("multicore", workers = 4)
options(future.globals.maxSize = 60 * 1024^3) # 60GB


# --------------------------------------------------------------------------
# Load the preprocess data
# --------------------------------------------------------------------------
seurat_file <- file.path(opt$input_dir, "preprocess_seurat.rds")
if (!file.exists(seurat_file)) stop("Error: Seurat object file not found: ", seurat_file)
seurat_object <- readRDS(seurat_file)
message("Seurat object loaded successfully.")

# --------------------------------------------------------------------------
# Function: Filter outlier cells
# --------------------------------------------------------------------------
filter_outliers <- function(sobj, lower_q, upper_q, max_mt, max_pt, log10GenesPerUMI_threshold) {
  lower_gene <- quantile(sobj@meta.data$nFeature_RNA, lower_q)
  lower_UMI <- quantile(sobj@meta.data$nCount_RNA, lower_q)
  upper_gene <- quantile(sobj@meta.data$nFeature_RNA, upper_q)
  
  sobj <- subset(sobj, subset =
                   (nCount_RNA >= lower_UMI) &
                   (nFeature_RNA >= lower_gene & nFeature_RNA <= upper_gene) &
                   (log10GenesPerUMI >= log10GenesPerUMI_threshold))
  
  if ("percent.mt" %in% colnames(sobj@meta.data)) {
    sobj <- subset(sobj, subset = percent.mt < max_mt)
  }
  if ("percent.pt" %in% colnames(sobj@meta.data)) {
    sobj <- subset(sobj, subset = percent.pt < max_pt)
  }
  
  return(sobj)
}

# --------------------------------------------------------------------------
# Apply filtering to each sample
# --------------------------------------------------------------------------
seurat_list <- SplitObject(seurat_object, split.by = "Sample")
seurat_list <- lapply(seurat_list, filter_outliers, lower_q = opt$lower_quantile,
                      upper_q = opt$upper_quantile, max_mt = opt$max_mt_percent,
                      max_pt = opt$max_pt_percent, log10GenesPerUMI_threshold=opt$log10GenesPerUMI_threshold)

# --------------------------------------------------------------------------
# Normalize data
# --------------------------------------------------------------------------
if (opt$normalization_method == "SCT") {
  message("Performing SCTransform normalization...")
  seurat_list <- lapply(X = seurat_list, FUN = SCTransform, method = "glmGamPoi")
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = opt$nfeatures, 
                                        assay = rep("SCT", length(seurat_list)))
  seurat_object <- merge(seurat_list[[1]], y = seurat_list[-1])
} else { # LogNormalize
  message("Performing LogNormalize...")
  seurat_list <- lapply(seurat_list, function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, nfeatures = opt$nfeatures)
    return(x)
  })
  seurat_object <- merge(seurat_list[[1]], y = seurat_list[-1])
  seurat_object <- FindVariableFeatures(seurat_object, nfeatures = opt$nfeatures)
  features <- VariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = features)
}

# --------------------------------------------------------------------------
# Merge samples and perform PCA
# --------------------------------------------------------------------------
message("Running PCA...")
seurat_object <- RunPCA(
  seurat_object,
  features = features,
  npcs = 100,
  verbose = FALSE
)

# --------------------------------------------------------------------------
# Save results
# --------------------------------------------------------------------------
ggsave(filename = file.path(output_dir, "elbow_plot_50.pdf"), plot = ElbowPlot(seurat_object, ndims = 50))
ggsave(filename = file.path(output_dir, "elbow_plot_100.pdf"), plot = ElbowPlot(seurat_object, ndims = 100))

saveRDS(seurat_object, file.path(output_dir, "filter_normailze_seurat.rds"))
message("filter normaization and pca completed successfully! Results saved in: ", output_dir) 
