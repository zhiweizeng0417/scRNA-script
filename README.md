# Single-cell RNA Sequencing Analysis Pipeline

This repository contains a suite of R and Python scripts designed for a comprehensive single-cell RNA sequencing (scRNA-seq) data analysis pipeline. The pipeline covers essential steps from raw data preprocessing to advanced analyses like cell type annotation and pseudotime trajectory inference.

## Scripts Overview

Here's a brief description of each script included in this repository:

* `scrublet_mult.py`: A Python script for doublet detection using Scrublet, processing multiple single-cell matrices.
* `scRNAseq_preprocessing.R`: R script for initial scRNA-seq data preprocessing, including loading, quality control metrics calculation (e.g., mitochondrial and ribosomal gene percentages), and basic filtering.
* `scRNAseq_filtering_normalization_pca.R`: R script to perform advanced filtering based on gene and UMI counts, normalize data using methods like SCTransform or LogNormalize, and execute Principal Component Analysis (PCA).
* `scRNAseq_Integrated.R`: R script for integrating multiple scRNA-seq samples, primarily using Seurat's `IntegrateData` function with SCTransform normalization, followed by PCA, UMAP, and t-SNE for integrated visualization.
* `scRNAseq_clustering_reduction.R`: R script dedicated to clustering and further dimensionality reduction (UMAP, t-SNE). It includes options for batch effect correction using Harmony or Seurat's integration.
* `scRNAseq_get_cluster_markers.R`: R script for identifying marker genes for each cell cluster. It calculates average expression and generates heatmaps of top marker genes.
* `scRNAseq_annotate_celltypes.R`: R script for annotating cell types based on identified marker genes. It generates visualizations like dot plots and heatmaps to support annotation.
* `scRNAseq_monocle2_step1.R`: R script for single-cell pseudotime trajectory analysis using Monocle2, including dimension reduction and cell ordering, and generating trajectory plots.

## Installation

### Prerequisites

* R (>= 4.0)
* Python (>= 3.7)
* Required R packages: `Seurat`, `tidyverse`, `patchwork`, `harmony`, `clustree`, `future`, `optparse`, `pheatmap`, `scRNAtoolVis`, `readxl`, `Biostrings`, `monocle`
* Required Python packages: `scrublet`, `numpy`, `scipy`, `matplotlib`, `pandas`, `numba`

### Setting up the Environment

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/zhiweizeng0417/scRNA-script
    cd scRNA-script
    ```

2.  **Install R packages:**
    Open R and run:
    ```R
    install.packages(c("tidyverse", "patchwork", "optparse", "pheatmap", "readxl", "Biostrings", "future",Seurat", "harmony","clustree", "monocle"))
    #For scRNAtoolVis, you might need to install from GitHub:https://github.com/junjunlab/scRNAtoolVis
    devtools::install_github("sajuukLyu/ggunchull", type = "source")
    devtools::install_github('junjunlab/scRNAtoolVis')
    ```

3.  **Install Python packages:**
    ```bash
    pip install scrublet numpy scipy matplotlib pandas numba
    ```

## Usage

Each script is designed to be run from the command line using `optparse` for R scripts and `argparse` for Python scripts. You can find detailed argument descriptions within each script by running them with the `--help` flag (e.g., `Rscript scRNAseq_preprocessing.R --help` or `python scrublet_mult.py --help`).

A typical workflow might look like this:

1.  **Doublet Detection:**
    ```bash
    python scrublet_mult.py -s /path/to/matrix_list.txt -o ./scrublet_results
    ```
     `matrix_list.txt` contains paths to your individual sample matrix directories：\
     /public/home/off_zengzhiwei/scRNA/LAp/matrix_data/YL \
     /public/home/off_zengzhiwei/scRNA/LAp/matrix_data/OL \
     /public/home/off_zengzhiwei/scRNA/LAp/matrix_data/YS \
     /public/home/off_zengzhiwei/scRNA/LAp/matrix_data/OS
 
3.  **Preprocessing:**
    ```bash
    Rscript scRNAseq_preprocessing.R --sample_names Sample1,Sample2 --matrix_dir /path/to/raw_data --doublet_dir ./scrublet_results --output_dir ./preprocess_results
    ```

4.  **Filtering, Normalization, PCA:**
    ```bash
    Rscript scRNAseq_filtering_normalization_pca.R --input_dir ./preprocess_results --output_dir ./norm_pca_results --normalization_method SCT --max_mt_percent 15
    ```
    
5.  **Clustering and Reduction:**
    ```bash
    Rscript scRNAseq_clustering_reduction.R --input_dir ./norm_pca_results --output_dir ./cluster_reduction_results --batch_correction harmony --pc_num 30 --cluster_resolutions 0.5,0.8,1.0
    ```

6.  **Get Cluster Markers:**
    ```bash
    Rscript scRNAseq_get_cluster_markers.R --input_dir ./cluster_reduction_results --output_dir ./marker_results --resolution 0.8 --top_n_markers 10
    ```

7.  **Cell Type Annotation:**
    ```bash
    Rscript scRNAseq_annotate_celltypes.R --input_dir ./marker_results --output_dir ./celltype_results --marker_genes CD3G,CD14,MS4A1 --cluster_order 5,3,1,2,4
    ```

8.  **Monocle2 Pseudotime Analysis:**
    ```bash
    Rscript scRNAseq_monocle2_step1.R --input_file /path/to/seurat_object.rds --output_dir ./monocle2_results --celltype_select "T cells" --root_celltype "Naive T cells"
    ```

Remember to adjust input and output directories, and other parameters, according to your specific data and analysis goals.
