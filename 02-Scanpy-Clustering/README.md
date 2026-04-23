# Part 2: Basic scRNA-seq Analysis — Preprocessing & Clustering

## Overview
This section performs a complete single-cell RNA-seq analysis workflow on the PBMC count matrix produced in Part 1, using **Scanpy** in Python.

## Input
Output files from Part 1:
- `../01-10X-Preprocessing/results/barcodes.tsv`
- `../01-10X-Preprocessing/results/genes.tsv`
- `../01-10X-Preprocessing/results/matrix.mtx`

## Notebook
`part2_scanpy_clustering.ipynb` — 16-step complete analysis notebook

## Pipeline Steps
| Step | Method | Result |
|------|--------|--------|
| Data Loading | scipy.io + AnnData | 252 cells × 2392 genes |
| Quality Control | calculate_qc_metrics | Ribosomal % flagged |
| Cell/Gene Filtering | filter_cells, filter_genes | Low quality removed |
| Doublet Detection | Scrublet | Doublets flagged |
| Normalization | normalize_total + log1p | Variance stabilized |
| Feature Selection | highly_variable_genes | 500 HVGs selected |
| PCA | tl.pca | 50 components |
| UMAP | tl.umap | 2D visualization |
| Clustering | Leiden algorithm | Multiple resolutions |
| Marker Genes | rank_genes_groups (Wilcoxon) | DEGs per cluster |
| Cell Annotation | Manual | PBMC cell types |

## How to Run
1. Open `part2_scanpy_clustering.ipynb` in Google Colab or Jupyter
2. Run Step 0 to install dependencies
3. Run Step 2 to upload the 3 MTX files from Part 1
4. Run all remaining cells sequentially
5. Download `pbmc_chrX_analyzed.h5ad` at the end

## Output Files (in `results/`)
| File | Description |
|------|-------------|
| `pbmc_chrX_analyzed.h5ad` | Fully analyzed AnnData object → input for Part 3 |
| `pbmc_cell_annotations.csv` | All cell metadata as CSV table |

## Reference
- [scverse Preprocessing and Clustering Tutorial](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/clustering.html)
