# --------------------------------------------------------------------------
# Script Name: scRNAseq_Integrated.R
# Author: Zeng Zhiwei
# Date: 2025-04-20
# Description: This script performs Integrated for single-cell RNA-seq data.
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
  make_option(c("-o", "--output_dir"), type = "character", default = "scRNAseq_Preprocess",
              help = "Output directory for results."),
  make_option(c("--celltype_file") ,type= "character" ,default = NULL,
              help = "celltype annotation information for per barcodes")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
# --------------------------------------------------------------------------
validate_parameters <- function(opt) {
  # Validate that required arguments are provided
  required_args <- c("sample_names", "matrix_dir")
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
}

validate_parameters(opt)

# Define paths and parameters
sample_names <- unlist(strsplit(opt$sample_names, ","))
matrix_dir <- opt$matrix_dir
doublet_dir <- opt$doublet_dir
output_dir <- opt$output_dir
celltype_dir <- opt$celltype_file
celltype_file<- read.csv(celltype_dir,header = T)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------
# Process each sample
# --------------------------------------------------------------------------
seurat_list <- list()
metadata <- c()

for (sample in sample_names) {
  sample_matrix_path <- file.path(matrix_dir, sample)
  
  if (!dir.exists(sample_matrix_path)) stop(paste("Matrix directory missing:", sample_matrix_path))
  
  seurat_list[[sample]] <- Read10X(data.dir = sample_matrix_path)
  
  total_cells <- ncol(seurat_list[[sample]])
  
  metadata <- c(metadata, rep(sample, total_cells))
  
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
# filter cell and intergated
# --------------------------------------------------------------------------
seurat_object <- subset(seurat_object,cells=celltype_file$barcode)
seurat_object$celltype<-celltype_file$celltype
# perform SCTransform normalization
seurat_list <- SplitObject(seurat_object, split.by = "Sample")
seurat_list <- lapply(X = seurat.list, FUN = SCTransform)
# select integration features and prep step
features <- SelectIntegrationFeatures(seurat_list)
seurat_list <- PrepSCTIntegration(
  seurat_list,
  anchor.features = features
)

# downstream integration steps
anchors <- FindIntegrationAnchors(
  seurat_list,
  normalization.method = "SCT",
  anchor.features = features
)
seurat.integrated <- IntegrateData(anchors, normalization.method = "SCT")

# PCA and UMAP 
seurat_object <- RunPCA(seurat.integrated, npcs = 30, verbose = FALSE)
seurat_object <- RunUMAP(seurat.integrated, reduction = 'pca', dims = 1:30)
seurat_object <- RunTSNE(seurat.integrated, reduction = 'pca', dims = 1:30)


p_umap_celltype <- DimPlot(seurat.integrated, reduction = "umap", group.by = 'celltype')
p_umap_sample <- DimPlot(seurat.integrated, reduction = "umap", group.by = 'Sample')
p_tsne_celltype <- DimPlot(seurat.integrated, reduction = "tsne", group.by = 'celltype')
p_tsne_sample <- DimPlot(seurat.integrated, reduction = "tsne", group.by = 'Sample')


ggsave(file.path(output_dir,"Integrated_umap_plot_sample.pdf"), plot = p_umap_sample)
ggsave(file.path(output_dir,"Integrated_umap_plot_celltype.pdf"), plot = p_umap_celltype)

ggsave(file.path(output_dir,"Integrated_tsne_plot_sample.pdf"), plot = p_tsne_sample)
ggsave(file.path(output_dir,"Integrated_tsne_plot_celltype.pdf"), plot = p_tsne_celltype)



# --------------------------------------------------------------------------
# Save the Seurat object
# --------------------------------------------------------------------------
saveRDS(seurat_object, file = file.path(output_dir, 'Integrated_seurat.rds'))
print("Integrated completed successfully.")
