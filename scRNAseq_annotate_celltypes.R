# --------------------------------------------------------------------------
# Script Name: annotate_celltypes.R
# Author: Zeng Zhiwei
# Date: 2025-03-23
# Description: This script loads a Seurat object with clustering information,
#              performs cell type annotation, generates visualizations, and
#              saves the annotated Seurat object.
# --------------------------------------------------------------------------

# Load required R packages
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(scRNAtoolVis)
  library(optparse)
  library(pheatmap)
})

# --------------------------------------------------------------------------
# Set up command line arguments
# --------------------------------------------------------------------------
option_list <- list(
  make_option("--input_dir", type = "character", default = NULL,
              help = "Directory containing the Seurat object with clustering information (required)."),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Output directory for annotation results (optional)."),
  make_option("--marker_genes", type = "character", default = NULL,
              help = "Comma-separated list of marker genes for dotplot."),
  make_option("--cluster_order", type = "character", default = NULL,
              help = "Comma-separated list of new cluster order."),
  make_option("--cell_types", type = "character", default = NULL,
              help = "Comma-separated list of new cell type identities."),
  make_option("--umap_colors", type = "character", default = NULL,
              help = "Comma-separated list of colors for UMAP plot.")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required arguments
if (is.null(opt$input_dir)) stop("Error: Input directory (--input_dir) is required.")

# Set up directories and parameters
input_dir <- normalizePath(opt$input_dir, mustWork = TRUE)
marker_genes_input <- opt$marker_genes
cluster_order_input <- opt$cluster_order
cell_types_input <- opt$cell_types
umap_colors_input <- opt$umap_colors

marker_genes <- if (!is.null(marker_genes_input)) strsplit(marker_genes_input, ",")[[1]] else NULL
cluster_order <- if (!is.null(cluster_order_input)) as.numeric(strsplit(cluster_order_input, ",")[[1]]) else NULL
cell_types <- if (!is.null(cell_types_input)) strsplit(cell_types_input, ",")[[1]] else NULL
umap_colors <- if (!is.null(umap_colors_input)) strsplit(umap_colors_input, ",")[[1]] else c(
  "#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2",
             "#8491B3B2", "#91D1C2B2", "#DC0000B2", "#7E61D0B2", "#374E55B2",
             "#DF8F44B2", "#A20056B2", "#6A6599B2", "#CB9A66B2", "#E69F00B2",
             "#56B4E9B2", "#800080B2", "#228B22B2", "#FFD700B2", "#008000B2",
             "#FF69B4B2", "#4682B4B2", "#FFA07AB2", "#9400D3B2", "#00CED1B2",
             "#F08080B2", "#00FF7FB2", "#DDA0DDB2", "#40E0D0B2", "#FF4500B2"
)

# Validate the length of cell types and cluster order if both are provided
if (!is.null(cell_types) && !is.null(cluster_order) && length(cell_types) != length(cluster_order)) {
  stop("Error: The number of cell type identities must match the number of clusters in the new order.")
}

# Define output directory
output_dir <- ifelse(is.null(opt$output_dir),
                     file.path(format(Sys.Date(), "%Y%m%d"), "celltype_annotation"),
                     normalizePath(opt$output_dir, mustWork = FALSE))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load Seurat object
seurat_file <- file.path(input_dir, "Seurat_0.4.rds")
if (!file.exists(seurat_file)) stop("Error: Seurat object not found: ", seurat_file)
seurat_object <- readRDS(seurat_file)
message("Seurat object loaded successfully.")

# --------------------------------------------------------------------------
# Reorder clusters (using input parameter)
# --------------------------------------------------------------------------
if (!is.null(cluster_order)) {
  seurat_object$seurat_clusters <- factor(seurat_object$seurat_clusters, levels = cluster_order)
  Idents(seurat_object) <- "seurat_clusters" # Make sure identities are based on the reordered clusters
  message("Clusters reordered based on provided parameter.")
} else {
  message("Note: Cluster order not provided, skipping cluster reordering.")
}

# --------------------------------------------------------------------------
# Plot marker genes (dotplot) (using input parameter)
# --------------------------------------------------------------------------
if (!is.null(marker_genes)) {
  p <- jjDotPlot(object =seurat_object , xtree = FALSE, ytree = FALSE, gene = marker_genes,
                 rescale = FALSE, dot.col = c('blue', 'white', 'red'), assay = 'SCT', midpoint = 0)
  ggsave(filename = file.path(output_dir, paste0('marker_dotplopt.pdf')), plot = p)
  message("Marker gene dotplot saved successfully in: ", output_dir)
} else {
  message("Note: Marker genes not provided, skipping dotplot.")
}

# --------------------------------------------------------------------------
# Set new identities (cell type annotation) (using input parameter)
# --------------------------------------------------------------------------
if (!is.null(cell_types) && !is.null(cluster_order)) {
  new_idents <- cell_types
  names(new_idents) <- levels(seurat_object) # Ensure names match the current cluster levels
  
  project <- RenameIdents(seurat_object, new_idents)
  project$celltype <- Idents(project)
  project$celltype <- factor(project$celltype,
                             levels = unique(new_idents)) # Use unique new_idents for levels
  message("Cell types annotated based on provided parameter.")
} else {
  message("Note: Cell type identities or cluster order not provided, skipping cell type annotation.")
  project <- seurat_object # Assign original object if annotation is skipped
}

# --------------------------------------------------------------------------
# Output UMAP plot with cell type labels (using input parameter)
# --------------------------------------------------------------------------
 

  
umap_plot_celltype <- clusterCornerAxes(object = project, reduction = 'umap',
                                        noSplit = TRUE,
                                        cornerTextSize = 3.5,
                                        addCircle = TRUE,
                                        cicAlpha = 0.1,
                                        nbin = 500,
                                        cicDelta = 0.2, cicLineSize = NA,
                                        pSize = 0.1, clusterCol = 'celltype',
                                        cellLabel = FALSE) +
  scale_color_manual(values = umap_colors) +
  scale_fill_manual(values = umap_colors)
ggsave(filename = file.path(output_dir, paste0("_umap_celltype.pdf")),
       plot = umap_plot_celltype, width = 8, height = 6)
message("UMAP plot with cell type labels saved successfully in: ", output_dir)

# --------------------------------------------------------------------------
# Plot cell type proportions
# --------------------------------------------------------------------------
cluster_freq <- as.data.frame(table(project@meta.data$celltype,
                                    project@meta.data$Sample))
cluster_freq <- cluster_freq %>%
  dplyr::group_by(Var2) %>%
  dplyr::mutate('percent' = Freq / sum(Freq))
colnames(cluster_freq) <- c('celltype', 'sample', 'count', 'percent')
p1 <- ggplot(cluster_freq, aes(x = celltype, y = percent, group = sample, color = sample)) +
  geom_point(size = 3) +
  geom_line(size = 1.5) +
  ylab(label = 'percent of celltype') +
  theme_minimal() +
  scale_color_manual(values = c('#6A5ACD', '#9370DB', '#33A02C')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
ggsave(filename = file.path(output_dir, paste0("celltype_percent.pdf")),
       plot = p1, width = 8, height = 5)

cluster_freq2 <- as.data.frame(table(seurat_object@meta.data$seurat_clusters,
                                     seurat_object@meta.data$Sample))
cluster_freq2 <- cluster_freq2 %>%
  dplyr::group_by(Var2) %>%
  dplyr::mutate('percent' = Freq / sum(Freq))
colnames(cluster_freq2) <- c('cluster', 'sample', 'count', 'percent')
p2 <- ggplot(cluster_freq2, aes(x = cluster, y = percent, group = sample, color = sample)) +
  geom_point(size = 3) +
  geom_line(size = 1.5) +
  ylab(label = 'percent of celltype') +
  theme_minimal() +
  scale_color_manual(values = c('#6A5ACD', '#9370DB', '#33A02C')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
ggsave(filename = file.path(output_dir, paste0("cluster_percent.pdf")),
       plot = p2, width = 10, height = 5)
message("Cell type proportion plots saved successfully in: ", output_dir)

# --------------------------------------------------------------------------
# Redraw heatmap after cell type annotation (assuming 'new_id' is available)
# --------------------------------------------------------------------------
if (!is.null(cell_types)) {
  exp_cluster_celltype <- AverageExpression(project, assays = 'SCT', group.by = 'celltype')
  cm_celltype <- exp_cluster_celltype$SCT
  Idents(seurat_object) <- Idents(project)
  seurat_object <- PrepSCTFindMarkers(seurat_object)
  findallmarker_celltype <- FindAllMarkers(seurat_object, logfc.threshold = 0.25,
                                           only.pos = TRUE,
                                           min.pct = 0.1, assay = "SCT", slot = 'data')
  
  # Filter genes
  diffs_celltype <- findallmarker_celltype %>%
    dplyr::filter(!grepl("^gene-Capa", gene, ignore.case = TRUE)) %>%
    dplyr::filter(avg_log2FC > 0.25, pct.1 > 0.25, pct.1 > pct.2) %>%
    distinct(gene, .keep_all = TRUE) %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(5, wt = avg_log2FC) %>%
    dplyr::arrange(cluster, -avg_log2FC)
  
  # Extract top 5 marker genes per cell type for heatmap
  top5_markers_celltype <- diffs_celltype %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 5, wt = avg_log2FC)
  
  cm_celltype <- cm_celltype[top5_markers_celltype$gene,]
  cm_celltype <- t(apply(cm_celltype, 1, scale))
  colnames(cm_celltype) <- levels(project)
  # Use the Gene_name_unique from top5_markers_celltype if available, otherwise use the original gene name
  row.names(cm_celltype) <- if ("Gene_name_unique" %in% colnames(top5_markers_celltype)) {
    top5_markers_celltype$Gene_name_unique
  } else {
    top5_markers_celltype$gene
  }
  
  pheatmap::pheatmap(cm_celltype, cluster_cols = FALSE, cluster_rows = FALSE,
                     color = colorRampPalette(c("#442e63", "#088682", "#f4f262"))(100),
                     fontsize_row = 7, border_color = NA, width = 5.5, height = 8, show_rownames = TRUE,
                     filename = file.path(output_dir, paste0('_heatmap_celltype.pdf')))
  message("Heatmap after cell type annotation saved successfully in: ", output_dir)
} else {
  message("Note: 'new_id' object not found. Skipping heatmap after cell type annotation.")
}


#cpm exp
exp_cluster <- AggregateExpression(project,group.by = c('celltype','Sample'),assays = 'SCT')
exp_cluster <- as.data.frame(exp_cluster$SCT)
total_counts<- colSums(exp_cluster)
cpm_exp_cluster <- t(apply(exp_cluster, 1, function(gene_counts) {
  (gene_counts / total_counts) * 1e6
}))

write.csv(cpm_exp_cluster, file = file.path(output_dir, paste0("CPM_exp_celltype.csv")), row.names = TRUE)

# --------------------------------------------------------------------------
# Save final Seurat object with cell type information
# --------------------------------------------------------------------------
saveRDS(project, file = file.path(output_dir, paste0('celltype.rds'))) # 保存 project 对象
message("Final Seurat object with cell type information saved successfully in: ", output_dir)