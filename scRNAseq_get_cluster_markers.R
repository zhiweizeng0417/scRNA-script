# --------------------------------------------------------------------------
# Script Name: get_cluster_markers.R
# Author: Zeng Zhiwei
# Date: 2025-03-22
# Description: This script loads clustered scRNA-seq data, identifies marker
#              genes for each cluster, calculates average expression, and
#              generates a top marker gene heatmap, UMAP, and t-SNE plots.
# --------------------------------------------------------------------------

# Load required R packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(pheatmap)
  library(optparse)
  library(scRNAtoolVis) # Make sure this package is installed
})

# --------------------------------------------------------------------------
# Set up command line arguments
# --------------------------------------------------------------------------
option_list <- list(
  make_option("--input_dir", type = "character", default = NULL,
              help = "Directory containing the clustered and reduced data (required)."),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Output directory for marker gene results (optional)."),
  make_option("--resolution", type = "numeric", default = 0.5,
              help = "Clustering resolution to use [default= %default]."),
  make_option("--top_n_markers", type = "integer", default = 5,
              help = "Number of top markers to display in the heatmap per cluster [default= %default].")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$input_dir)) stop("Error: Input directory (--input_dir) is required.")

# Set up directories and parameters
input_dir <- normalizePath(opt$input_dir, mustWork = TRUE)
resolution <- opt$resolution
top_n_markers <- opt$top_n_markers

# Define output directory
if (is.null(opt$output_dir)) {
  output_dir <- file.path("scRNAseq_results", "cluster_marker", format(Sys.Date(), "%Y%m%d"))
} else {
  output_dir <- normalizePath(opt$output_dir, mustWork = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load Seurat object
seurat_file <- file.path(input_dir, "clustered_reduced_data.rds")
if (!file.exists(seurat_file)) stop("Error: Clustered and reduced data file not found: ", seurat_file)
seurat_object <- readRDS(seurat_file)
message("Clustered and reduced data loaded successfully.")

# --------------------------------------------------------------------------
# Set resolution and update cluster identities
# --------------------------------------------------------------------------
cluster_column <- paste0("SCT_snn_res.", resolution)
if (!cluster_column %in% colnames(seurat_object@meta.data)) {
  stop(paste("Error: Clustering resolution '", resolution, "' not found in Seurat object metadata."))
}
seurat_object@meta.data$seurat_clusters <- seurat_object@meta.data[[cluster_column]]
Idents(seurat_object) <- "seurat_clusters"
seurat_object$seurat_clusters <- as.numeric(as.character(seurat_object$seurat_clusters))
seurat_object$seurat_clusters <- factor(seurat_object$seurat_clusters, levels = sort(unique(seurat_object$seurat_clusters)))

# --------------------------------------------------------------------------
# Define cluster colors (at least 30 colors)
# --------------------------------------------------------------------------
cluster_colors <- c(
  "#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2",
  "#8491B3B2", "#91D1C2B2", "#DC0000B2", "#7E61D0B2", "#374E55B2",
  "#DF8F44B2", "#A20056B2", "#6A6599B2", "#CB9A66B2", "#E69F00B2",
  "#56B4E9B2", "#800080B2", "#228B22B2", "#FFD700B2", "#008000B2",
  "#FF69B4B2", "#4682B4B2", "#FFA07AB2", "#9400D3B2", "#00CED1B2",
  "#F08080B2", "#00FF7FB2", "#DDA0DDB2", "#40E0D0B2", "#FF4500B2"
)

# Ensure that the number of colors is sufficient
if (length(cluster_colors) < length(unique(seurat_object$seurat_clusters))) {
  warning(paste("Warning: Number of colors provided (", length(cluster_colors), ") is less than the number of clusters (", length(unique(seurat_object$seurat_clusters)), "). Colors will be recycled."))
}

# --------------------------------------------------------------------------
# Output UMAP and t-SNE plots
# --------------------------------------------------------------------------
# UMAP plot output
umap_plot <- clusterCornerAxes(object = seurat_object, reduction = 'umap', noSplit = TRUE, 
                               arrowType = 'open', 
                               pSize = 0.1) +
  scale_color_manual(values = cluster_colors)
ggsave(filename = file.path(output_dir, paste0("umap_cluster.pdf")), plot = umap_plot, width = 8, height = 6)
message("UMAP plot colored by cluster saved successfully in: ", output_dir)

# t-SNE plot output
tsne_plot <- clusterCornerAxes(object = seurat_object, reduction = 'tsne', noSplit = TRUE,
                               arrowType = 'open',
                               pSize = 0.1, 
                               cellLabel = FALSE) +
  scale_color_manual(values = cluster_colors)
ggsave(filename = file.path(output_dir, paste0("tsne_cluster.pdf")), plot = tsne_plot, width = 8, height = 6)
message("t-SNE plot colored by cluster saved successfully in: ", output_dir)

# --------------------------------------------------------------------------
# Calculate marker genes for each cluster
# --------------------------------------------------------------------------
# No need for PrepSCTFindMarkers anymore
seurat_object<-PrepSCTFindMarkers(seurat_object,assay = 'SCT')
findallmarker <- FindAllMarkers(seurat_object, logfc.threshold = 0.25, 
                                only.pos = TRUE, test.use = "bimod", min.pct = 0.1)

# Filter genes
diffs <- findallmarker %>%
  dplyr::filter(!grepl("^gene-Capa", gene, ignore.case = TRUE)) %>%
  dplyr::filter(avg_log2FC > 0.25, pct.1 > 0.25, pct.1 > pct.2) %>%
  distinct(gene, .keep_all = TRUE)

# Extract top 50 marker genes per cluster
top50 <- diffs %>%
  dplyr::group_by(cluster) %>%
  top_n(50, wt = avg_log2FC) %>%
  dplyr::arrange(cluster, -avg_log2FC)

#cpm exp
exp_cluster <- AggregateExpression(seurat_object,group.by = 'seurat_clusters',assays = 'SCT')
exp_cluster <- as.data.frame(expr_cluster$SCT)
total_counts<- colSums(exp_cluster)
cpm_exp_cluster <- t(apply(exp_cluster, 1, function(gene_counts) {
  (gene_counts / total_counts) * 1e6
}))


# Save marker gene tables
write.csv(diffs, file = file.path(output_dir, paste0('findallmarker.csv')), row.names = FALSE)
write.csv(cpm_exp_cluster, file = file.path(output_dir, paste0("CPM_exp_cluster.csv")), row.names = TRUE)
write.csv(top50, file = file.path(output_dir, paste0("top50marker.csv")), row.names = FALSE)
message("Marker gene tables saved successfully in: ", output_dir)

# --------------------------------------------------------------------------
# Generate top N heatmap
# --------------------------------------------------------------------------
# Use the default assay for the heatmap
cm <- cpm_exp_cluster   
top_markers_heatmap <- top50 %>%
  dplyr::group_by(cluster) %>%
  top_n(n = top_n_markers, wt = avg_log2FC) %>%
  dplyr::arrange(cluster, -avg_log2FC)
cm_heatmap <- cm[top_markers_heatmap$gene,]
cm_heatmap <- t(apply(cm_heatmap, 1, scale))
colnames(cm_heatmap) <- levels(seurat_object)
cm_heatmap[cm_heatmap > 3] <- 3
pheatmap::pheatmap(cm_heatmap, cluster_cols = FALSE, cluster_rows = FALSE, color = colorRampPalette(c("#442e63", "#088682", "#f4f262"))(100),
                   fontsize_row = 7, border_color = NA, width = 5.5, height = 8, show_rownames = TRUE,
                   filename = file.path(output_dir, paste0('heatmap_top', top_n_markers, '.pdf')))
message("Top ", top_n_markers, " marker gene heatmap saved successfully in: ", output_dir)

# Save the Seurat object with marker information (optional)
saveRDS(seurat_object, file.path(output_dir, paste0('Seurat_',resolution,'.rds')))
message("Seurat object with marker information saved (optional) in: ", output_dir)