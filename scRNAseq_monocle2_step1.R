# ------------------------------------------------------------------------------
# Script Name: monocle2 dimension reduce and ordering cell.R
# Author: Zeng Zhiwei
# Date: 2025-06-01
# Description: This script performs scRNAcell-seq pseudo-time trajectory analysis
#              with monocle2.
# ------------------------------------------------------------------------------

# Load required R packages
suppressPackageStartupMessages({
  library(Seurat,verbose = F)
  library(tidyverse,verbose = F)
  library(patchwork,verbose = F)
  library(monocle,verbose = F)
  library(future,verbose = F)
})

# ------------------------------------------------------------------------------
# Set up command line arguments
# ------------------------------------------------------------------------------
option_list <- list(
  make_option("--input_file", type = "character", default = NULL,
              help =  "Path to the Seurat object (.rds file) for cell trajectory 
              analysis."),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Output directory for cell trajectory analysis results."),
  make_option("--celltype_select", type = "character", default = NULL,
              help = "Cell type to performed cell trajectory analysis."),
  make_option("--root_celltype", type = "character", default = NULL,
              help = "Root cell type to sort cell"),
  make_option("--dim_use", type = "integer", default = 2,
              help = "Number of dim to use to perform DDRTree reduce."),
  make_option("--ordering_genes_file", type = "character", default = NULL,
              help = "Optional: Path to a plain text file (one gene per line) 
              containing a custom list of genes to use for cell ordering in 
              Monocle2. If not provided, differentially expressed genes 
              (via Seurat's FindAllMarkers) will be used."),
  make_option("--n_cores", type = "integer", default = 4,
              help = "Number of cores to use for parallel processing.")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# ------------------------------------------------------------------------------
# Argument Validation and Preparation
# ------------------------------------------------------------------------------
# Check for mandatory arguments if they are truly required
if (is.null(opt$input_file)) {
  stop("Error: --input_file is a mandatory argument. Please specify the input directory.")
}
if (is.null(opt$output_dir)) {
  stop("Error: --output_dir is a mandatory argument. Please specify the output directory.")
}
if (is.null(opt$celltype_select)) {
  stop("Error: --celltype_select is a mandatory argument. Please specify the output directory.")
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
  message("Created output directory: ", opt$output_dir)
}


# Set up directories and parameters
input_file <- normalizePath(opt$input_file, mustWork = TRUE)
output_dir <- opt$output_dir
celltype_select_vector <- unlist(strsplit(opt$celltype_select, ","))
root_celltype <- opt$root_celltype
dim_use<-opt$dim_use
ordering_genes_file <- opt$ordering_genes_file
n_cores <- opt$n_cores

# Set up parallel processing (e.g., for Monocle2 steps that support it)
options(mc.cores = n_cores) # For monocle's internal mclapply, if used
plan(multisession, workers = n_cores) # For future-based parallelization

message(paste0("Analysis will use ", n_cores, " cores for parallel processing."))


# ------------------------------------------------------------------------------
# Load Seurat object and Prepare Monocle2 CellDataSet
# ------------------------------------------------------------------------------

message("Loading Seurat object from: ", input_file)
if (!file.exists(input_file)) {
  stop("Error: Processed Seurat object file not found at specified path: ", input_file, call. = FALSE)
}
project <- readRDS(input_file)
message("Seurat object loaded successfully.")

# Subset Seurat object to selected cell types *before* Monocle2 conversion
message("Subsetting Seurat object to include only selected cell types: ", paste(celltype_select_vector, collapse = ", "))
project_subset <- project[, project@meta.data$cell_type %in% celltype_select_vector]

if (ncol(project_subset) == 0) {
  stop("Error: No cells found for the selected cell types. Please check --celltype_select values.", call. = FALSE)
}
message("Number of cells selected for analysis: ", ncol(project_subset))

# Extract data for Monocle2 CellDataSet
cm <- GetAssayData(project_subset, slot = 'counts') # Ensure you have a 'counts' assay and slot
md <- project_subset@meta.data
gene_annotation <- data.frame(gene_id = rownames(cm),
                              gene_short_name = rownames(cm)) # Monocle2 often likes gene_short_name
rownames(gene_annotation) <- as.vector(gene_annotation$gene_id) # Ensure rownames match gene_id
fd <- new("AnnotatedDataFrame", data = gene_annotation)
pd <- new("AnnotatedDataFrame", data = md)

# Create CellDataSet object
monocl_data <- newCellDataSet(cm,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = negbinomial.size())
message("Monocle2 CellDataSet object created.")

# Estimate size factors and dispersions
monocl_data <- estimateSizeFactors(monocl_data)
monocl_data <- estimateDispersions(monocl_data)
message("Size factors and dispersions estimated.")

# ------------------------------------------------------------------------------
# Gene Selection for Ordering Cells
# ------------------------------------------------------------------------------
ordering_genes <- NULL

if (!is.null(ordering_genes_file)) {
  # Option 1: Use a custom gene list from a file
  if (!file.exists(ordering_genes_file)) {
    stop("Error: Custom ordering genes file not found: ", ordering_genes_file, 
         call. = FALSE)
  }
  custom_genes <- read.table(ordering_genes_file,header = T)
  # Filter to only genes present in the data
  ordering_genes <- intersect(custom_genes$ordering_genes, rownames(fData(monocl_data)))
  if (length(ordering_genes) == 0) {
    stop("Error: No genes from the provided --ordering_genes_file were found in 
         the dataset. Please check the file and gene names.", call. = FALSE)
  }
  message(paste0("Using ", length(ordering_genes), " custom genes from '", 
                 ordering_genes_file, "' for ordering filter."))
  
} else {
  # Option 2: Use differentially expressed genes (your original logic)
  message("No custom ordering genes file provided. Identifying differentially 
          expressed genes for ordering filter.")
  project_subset <- SetIdent(project_subset, value = "cell_type") 
  deg_genes_df <- FindAllMarkers(project_subset,
                                 logfc.threshold = 0.25,
                                 only.pos = TRUE,
                                 test.use = "bimod",
                                 min.pct = 0.1,
                                 assay = "SCT",
                                 slot = 'data')
  diffs2_genes <- deg_genes_df %>%
    filter(avg_log2FC > 0.25 & pct.1 > 0.3 & (pct.1 > pct.2)) %>%
    pull(gene) # Get the gene names as a vector
  
  if (length(diffs2_genes) == 0) {
    stop("Error: No differentially expressed genes passed the filtering criteria
         (avg_log2FC > 0.25, pct.1 > 0.3, pct.1 > pct.2). Please adjust filtering 
         or provide a custom gene list.", call. = FALSE)
  }
  ordering_genes <- diffs2_genes
  message(paste0("Found ", length(ordering_genes), " differentially expressed genes for ordering filter."))
}

# Apply the selected ordering genes to Monocle2 object
monocl_data <- setOrderingFilter(monocl_data, ordering_genes)
message("Monocle2 ordering filter set.")

# ------------------------------------------------------------------------------
# Dimensionality Reduction (DDRTree)
# ------------------------------------------------------------------------------
message(paste0("Performing DDRTree dimensionality reduction using ", dim_use, " dimensions."))
monocl_data <- reduceDimension(monocl_data,
                               max_components = dim_use,
                               reduction_method = 'DDRTree',
                               norm_method = "log",
                               pseudo_count_expression = 1,
                               #residualModelFormulaStr = '~sample', # if you want remove sample infect
                               verbose = TRUE)
message("DDRTree dimensionality reduction complete.")

# ------------------------------------------------------------------------------
# Cell Ordering
# ------------------------------------------------------------------------------
message("Ordering cells by pseudo-time.")
if (!is.null(root_celltype)) {
  #order cell
  monocl_data <- orderCells(monocl_data)
  
  #select root cell associate state
  state_celltype_counts <- table(pData(monocl_data)$State, pData(monocl_data)$cell_type)
  counts_for_root_celltype <- state_celltype_counts[, root_celltype]
  root_state_to_use <- as.numeric(names(which.max(counts_for_root_celltype)))
  if (length(root_state_to_use) == 0 || is.na(root_state_to_use)) {
    stop(paste0("Error: Could not determine a clear root state for cell type '", 
                root_celltype, "'. No cells of this type were found in any Monocle2 
                state after DDRTree."), call. = FALSE)
  }
  
  #reorder cell 
  monocl_data <- orderCells(monocl_data, root_state = root_state_to_use)
  message(paste0("Using Monocle2 State '", root_state_to_use, "' as the root, 
                 based on its enrichment for '", root_celltype, "' cells."))
} else {
  monocl_data <- orderCells(monocl_data)
  message("Monocle2 inferred the root cell automatically.")
}
message("Cell ordering complete. Pseudo-time values assigned.")


# ------------------------------------------------------------------------------
# Cell Trajectory Plots
# ------------------------------------------------------------------------------
message("Generating cell trajectory plots...")
col=c("#E31A1C", "#A6CEE3", "#1F78B4", "#FB9A99", "#FDBF6F", "#FF7F00", "#CAB2D6", 
               "#6A3D9A", "#FFFF99", "#B15928", "#8DD3C7", "#FFFFB3")

# Plot 1: Cell Trajectory colored by Cell Type
p1 <- plot_cell_trajectory(monocl_data,
                           color_by = 'celltype',
                           cell_size = 0.1,
                           show_branch_points = T) +
  scale_color_manual(values = col)
  theme(legend.position = "right") 

ggsave(filename = file.path(output_dir, 'cell_trajectory_by_celltype.pdf'), 
       plot = p1, width = 6, height = 3.5)
message("Saved cell_trajectory_by_celltype.pdf")


# Plot 2: Cell Trajectory colored by Pseudotime
p2 <- plot_cell_trajectory(monocl_data,
                           color_by = 'Pseudotime',
                           cell_size = 0.1,
                           show_branch_points = T) +
  scale_color_gradient(low = "#ADD8E6", high = "#00008B") +
  theme(legend.position = "right") 

ggsave(filename = file.path(output_dir, 'cell_trajectory_pseudotime.pdf'), 
       plot = p2, width = 6, height = 3.5)
message("Saved cell_trajectory_pseudotime.pdf")

# Save final monocle2 object
saveRDS(monocl_data, file.path(output_dir,"monocl_data.rds"))
message("monocle2 cell trajectory successfully. Results saved in: ", output_dir)

