# Part 3: Getting Started with AnnData

## Overview
This section explores the **AnnData format** in depth using the analyzed PBMC dataset from Part 2. AnnData is the standard data structure for single-cell analysis in Python (scverse ecosystem).

## Input
- `../02-Scanpy-Clustering/results/pbmc_chrX_analyzed.h5ad`

## Notebook
`part3_anndata_tutorial.ipynb` — 15-step complete tutorial notebook

## AnnData Components Explored
| Component | What it stores | Our data |
|-----------|---------------|----------|
| `adata.X` | Active data matrix | Log-normalized counts |
| `adata.obs` | Cell metadata | Barcodes, QC, clusters, cell types |
| `adata.var` | Gene metadata | Ensembl IDs, HVG flags |
| `adata.obsm` | Cell embeddings | X_pca (50D), X_umap (2D) |
| `adata.layers` | Data versions | Raw, log-norm, CPM |
| `adata.uns` | Unstructured metadata | Project info, params |
| `adata.obsp` | Cell-cell matrices | k-NN connectivities |

## Key Concepts Demonstrated
- Loading h5ad files with `anndata.read_h5ad()`
- Inspecting all AnnData slots
- Subsetting by index, name, boolean mask
- Views vs copies (memory efficiency)
- Adding new layers (CPM normalization)
- Manual plotting from obsm embeddings
- Visualizing the connectivity matrix
- Exporting to DataFrames and CSV
- Writing enriched AnnData to disk

## How to Run
1. Open `part3_anndata_tutorial.ipynb` in Google Colab or Jupyter
2. Upload `pbmc_chrX_analyzed.h5ad` from Part 2 when prompted
3. Run all cells sequentially

## Output Files (in `results/`)
| File | Description |
|------|-------------|
| `pbmc_chrX_final.h5ad` | Final enriched AnnData with all pipeline results |

## References
- [Getting started with anndata](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
- [scverse AnnData getting started](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)
