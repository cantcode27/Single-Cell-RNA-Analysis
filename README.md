# Single-Cell RNA-Seq Analysis Pipeline

This repo covers a complete scRNA-seq analysis workflow broken into three self-contained notebooks — from raw 10X Genomics output all the way to clustered, annotated cell populations. Done as part of a Bioinformatics course assignment, using Python and Scanpy throughout.

The three modules build on each other sequentially, so if you're following along it makes sense to go in order.

---

## What's single-cell RNA-seq and why does it matter?

Regular bulk RNA-seq tells you the average gene expression across millions of cells in a sample — which is useful, but hides a lot. A tumor isn't one thing, it's thousands of different cell types doing different things. scRNA-seq lets you measure gene expression in **individual cells**, which means you can actually identify distinct cell populations, find rare cell types, and understand how cells differ from each other within the same tissue.

The tradeoff is complexity. The data is sparse (most genes aren't expressed in any given cell), noisy, and the preprocessing steps matter a lot. That's what this pipeline works through.

---

## Repo Structure

```
Single-Cell-RNA-Analysis/
│
├── 01-10X-Preprocessing/       # Loading and QC of raw 10X data
│   └── preprocessing_10x.ipynb
│
├── 02-Scanpy-Clustering/       # Dimensionality reduction and clustering
│   └── scanpy_clustering.ipynb
│
├── 03-AnnData-Tutorial/        # Deep dive into the AnnData object structure
│   └── anndata_tutorial.ipynb
│
└── README.md
```

---

## Module Breakdown

### 01 — 10X Genomics Preprocessing

10X Genomics Chromium is the most widely used scRNA-seq platform right now. After sequencing, Cell Ranger processes the raw reads and spits out three files: a barcode list, a gene list, and a sparse count matrix. This notebook picks up from there.

**What this notebook covers:**

- Loading the Cell Ranger output (the `filtered_feature_bc_matrix/` folder or equivalent) into a Scanpy AnnData object
- Basic QC metrics per cell: number of genes detected (`n_genes_by_counts`), total counts (`total_counts`), and mitochondrial gene percentage (`pct_counts_mt`)
- Filtering out low-quality cells — cells with too few genes are likely empty droplets, cells with too many genes might be doublets, and cells with high mitochondrial content are usually dying
- Visualizing QC distributions with violin plots and scatter plots to set reasonable thresholds
- Normalization: library-size normalization to 10,000 counts per cell followed by log1p transformation
- Identifying highly variable genes (HVGs) — the ~2,000 genes that vary most across cells, which are the ones that actually carry biological signal for clustering

By the end of this notebook you have a clean, normalized AnnData object ready for downstream analysis.

---

### 02 — Scanpy Clustering

This is where the biology starts coming through. Once you have a clean count matrix, the goal is to group cells that look similar (in terms of gene expression) together — and then figure out what those groups actually are biologically.

**What this notebook covers:**

- Scaling the data (zero mean, unit variance per gene) before PCA
- Running PCA and using an elbow plot / variance explained to decide how many principal components to keep
- Building a neighborhood graph on the PCA-reduced data — this is the graph that all downstream steps depend on
- UMAP embedding for visualization — projects the high-dimensional data down to 2D so you can actually look at it
- Leiden clustering (or Louvain, depending on the tutorial version) — community detection on the neighborhood graph that assigns each cell to a cluster
- Finding marker genes for each cluster using `sc.tl.rank_genes_groups()` — these are the genes that are significantly more expressed in one cluster vs. all others, which is how you figure out what cell type each cluster is
- Dot plots, heatmaps, and UMAP plots colored by cluster and by marker gene expression

The main output is a UMAP where each dot is a cell, colored by cluster — and a list of marker genes per cluster that you use to annotate what those clusters actually represent (e.g., "cluster 3 expresses CD3E and CD8A, so it's cytotoxic T cells").

---

### 03 — AnnData Tutorial

AnnData is the core data structure that Scanpy (and most of the Python single-cell ecosystem) is built around. If you don't understand it, you spend a lot of time confused about where things are stored and why operations work the way they do. This notebook is a focused tutorial on the object itself.

**What this notebook covers:**

- The AnnData structure: `adata.X` (the count matrix), `adata.obs` (cell-level metadata), `adata.var` (gene-level metadata), `adata.uns` (unstructured data like color palettes and PCA results), `adata.obsm` (multi-dimensional cell embeddings like UMAP coordinates), and `adata.layers` (storing multiple versions of the matrix — e.g., raw counts and normalized counts side by side)
- Creating AnnData objects from scratch
- Slicing and subsetting (rows = cells, columns = genes)
- Reading and writing `.h5ad` files (the standard format for saving AnnData objects)
- How Scanpy functions modify the AnnData in-place and where they store results

This one is more of a reference/tutorial than an analysis — but it's genuinely useful because once you understand the data structure everything else in Scanpy makes more sense.

---

## Tools & Dependencies

Everything runs in Python. Install dependencies with:

```bash
pip install scanpy anndata matplotlib seaborn pandas numpy
```

Or if you're using conda:

```bash
conda install -c conda-forge scanpy python-igraph leidenalg
```

Note: Leiden clustering requires `python-igraph` and `leidenalg` to be installed separately — they're not pulled in automatically by Scanpy.

| Library | Purpose |
|---------|---------|
| `scanpy` | Core scRNA-seq analysis toolkit |
| `anndata` | AnnData data structure |
| `numpy` / `pandas` | Data manipulation |
| `matplotlib` / `seaborn` | Plotting |
| `leidenalg` | Leiden clustering algorithm |
| `python-igraph` | Graph operations for clustering |

---

## How to Run

1. Clone the repo:
   ```bash
   git clone https://github.com/cantcode27/Single-Cell-RNA-Analysis.git
   cd Single-Cell-RNA-Analysis
   ```

2. Install dependencies (see above)

3. Open and run the notebooks in order:
   ```bash
   jupyter notebook
   ```
   Start with `01-10X-Preprocessing`, then `02-Scanpy-Clustering`, then `03-AnnData-Tutorial`.

The notebooks use publicly available demo datasets (PBMC 3k or similar) so you don't need to download anything separately — dataset loading is handled inside each notebook.

---

## The Broader Pipeline

```
Cell Ranger output (barcodes + genes + matrix)
         │
         ▼
  01 — Load → QC filtering → Normalize → Find HVGs
         │
         ▼
  02 — Scale → PCA → Neighbors → UMAP → Cluster → Marker genes
         │
         ▼
  03 — Understand the AnnData object that stores all of the above
```

---

## References

- Wolf FA, Angerer P, Theis FJ (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 19:15.
- Virshup I, et al. (2021). The scverse project provides a computational ecosystem for single-cell omics data analysis. *Nature Biotechnology*.
- 10X Genomics Cell Ranger documentation: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
- Scanpy tutorials: https://scanpy-tutorials.readthedocs.io/

---

*Bioinformatics course assignment — NUST SINES, Spring 2026.*
