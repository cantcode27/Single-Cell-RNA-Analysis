# Single-Cell RNA-Seq Analysis Pipeline 

This repo walks through a complete single-cell RNA-sequencing analysis pipeline on **human peripheral blood mononuclear cells (PBMCs)**, starting from raw FASTQ files and ending with annotated cell types stored in a final AnnData object. The work is split across three parts, each building directly on the previous one.

The dataset is a subset of 10X Genomics' 1k PBMC v3 library, aligned to **chromosome X only** (hg19), which keeps the gene space manageable while covering real immune cell biology.

---

## What's the biological question?

PBMCs are the immune cells circulating in blood — T cells, B cells, NK cells, monocytes, and so on. Single-cell RNA-seq lets you measure gene expression in each cell individually, so instead of seeing an average across everything, you can actually distinguish the different cell types and understand what each one is doing.

This pipeline takes that concept from raw sequencing reads all the way to a UMAP plot where each dot is a cell, colored by what type it is. The chrX-only constraint means we're working with 2,392 genes instead of ~30,000, but the methods are identical to a full genome analysis — it just runs faster.

---

## Repo Structure

```
Single-Cell-RNA-Analysis/
│
├── 01-10X-Preprocessing/
│   └── results/
│       ├── barcodes.tsv       
│       ├── genes.tsv          
│       └── matrix.mtx         
│
├── 02-Scanpy-Clustering/
│   └── results/
│       ├── pbmc_chrX_analyzed.h5ad     
│       └── pbmc_cell_annotations.csv    
│
├── 03-AnnData-Tutorial/
│   └── results/
│       └── pbmc_chrX_final.h5ad       
│
└── README.md
```

---

## Part 1 — 10X Genomics Preprocessing (Galaxy)

### What this part does

Raw FASTQ files from a 10X Chromium sequencing run don't come with a cell-by-gene count matrix — you have to build that yourself. This part handles that using two Galaxy tools: **RNA STARsolo** (alignment + quantification) and **DropletUtils** (cell calling).

This part runs entirely on [usegalaxy.org](https://usegalaxy.org) — no local installation needed. The input files come from [Zenodo record 3457880](https://zenodo.org/record/3457880).

### Input data

| File | What it is |
|------|------------|
| `subset_pbmc_1k_v3_S1_L001_R1_001.fastq.gz` | Lane 1 — cell barcodes + UMIs |
| `subset_pbmc_1k_v3_S1_L001_R2_001.fastq.gz` | Lane 1 — cDNA reads |
| `subset_pbmc_1k_v3_S1_L002_R1_001.fastq.gz` | Lane 2 — cell barcodes + UMIs |
| `subset_pbmc_1k_v3_S1_L002_R2_001.fastq.gz` | Lane 2 — cDNA reads |
| `3M-february-2018.txt.gz` | 10X v3 barcode whitelist (3 million barcodes) |
| `Homo_sapiens.GRCh37.75.gtf` | Gene annotation for hg19 |

### Pipeline steps

**Step 1 — RNA STARsolo**

STARsolo does two things simultaneously: it aligns the cDNA reads (R2) to the hg19 chrX reference genome, and it demultiplexes individual cells using the barcode whitelist (from R1). The result is a raw count matrix where every cell barcode gets a column and every gene gets a row. At this stage, 5,200 barcodes were detected — but most of these are empty droplets, not real cells.

**Step 2 — MultiQC**

A quick quality check on the STARsolo alignment log to confirm mapping rates look sensible before proceeding.

**Step 3 — DropletUtils (EmptyDrops)**

This is the key filtering step. In 10X sequencing, many droplets don't capture a cell at all — they just contain ambient RNA. DropletUtils uses a statistical test (EmptyDrops) to distinguish real cells from empty droplets, using a lower threshold of 200 UMIs and FDR of 0.01.

**Result: 5,200 raw barcodes → 252 high-quality cells.**

### Output files (in `results/`)

| File | Description |
|------|-------------|
| `barcodes.tsv` | 252 16-mer cell barcodes (one per line) |
| `genes.tsv` | 2,392 genes — Ensembl ID + gene name (tab-separated) |
| `matrix.mtx` | Sparse count matrix in MatrixMarket format (2392 genes × 252 cells, 39,002 non-zero entries) |

### Key numbers

| Metric | Value |
|--------|-------|
| Raw barcodes from STARsolo | 5,200 |
| High-quality cells after EmptyDrops | **252** |
| Genes quantified (chrX only) | 2,392 |
| Non-zero entries in matrix | 39,002 |

---

## Part 2 — Preprocessing & Clustering (Scanpy)

### What this part does

This is the core analysis notebook. It takes the 3 output files from Part 1 and runs a complete scRNA-seq preprocessing and clustering workflow in Python using **Scanpy**. By the end, every cell has a cluster assignment and a cell type label.

**Library versions used:** Scanpy 1.12.1 · AnnData 0.12.10 · NumPy 2.0.2 · Pandas 2.3.3

### Step-by-step walkthrough

**Step 3 — Load data into AnnData**

The three MTX files are loaded into an **AnnData** object — the standard data structure for single-cell analysis in Python. Genes from the `genes.tsv` file have two columns (Ensembl ID and gene name); gene names are used as the primary identifier.

```
AnnData object with n_obs × n_vars = 252 × 2392
```

**Step 4 — Quality Control**

Three core QC metrics are calculated per cell:
- `n_genes_by_counts` — number of unique genes detected per cell
- `total_counts` — total UMI count per cell
- `pct_counts_ribo` — percentage of reads mapping to ribosomal genes (54 ribosomal genes found on chrX)

QC summary across 252 cells:

| Metric | Min | Median | Max |
|--------|-----|--------|-----|
| Genes per cell | 78 | 146 | 451 |
| Total UMI counts | 245 | 1,039 | 7,780 |

Violin plots and scatter plots were used to inspect the distributions. Note: since this is a chrX-only dataset, there are no mitochondrial genes (MT- genes are on the mitochondrial chromosome, not chrX), so the usual `pct_counts_mt` filter isn't applicable here. Ribosomal content is flagged instead.

**Step 5 — Gene filtering**

Genes detected in fewer than 3 cells were removed. This dropped 1,641 genes, leaving **751 genes** that are expressed in at least 3 cells. Cell count stays at 252.

```
filtered out 1,641 genes detected in fewer than 3 cells
Cells: 252 | Genes: 751
```

**Step 5 — Doublet Detection (Scrublet)**

Doublets are droplets that accidentally captured two cells — they look like cells with unusually high gene counts and can distort clustering. Scrublet simulates synthetic doublets and gives each cell a `doublet_score` (0–1). Cells above the threshold are flagged with `predicted_doublet = True`. Doublets are flagged but kept in the dataset for Part 3's AnnData exploration; they're visible in the UMAP as outlier cells.

**Step 6 — Normalization**

Two-step normalization:
1. Raw counts saved to `adata.layers["counts"]` for safekeeping
2. Each cell's counts normalized to the median total count depth (library-size normalization)
3. Log1p transformation: `X_normalized = log(1 + X)` — stabilizes variance so highly expressed genes don't dominate downstream analysis

**Step 7 — Highly Variable Genes (HVGs)**

Not all 751 genes carry useful information for distinguishing cell types. HVG selection identifies the genes with the most variability across cells — these are the ones that actually differentiate one cell population from another.

```
Number of highly variable genes selected: 500
```

Top HVGs by normalized dispersion:

| Gene | Ensembl ID | Dispersion (norm) |
|------|------------|-------------------|
| TSPYL2 | ENSG00000184205 | 4.93 |
| IL3RA | ENSG00000185291 | 3.41 |
| SH2D1A | ENSG00000183918 | 3.25 |
| LINC00892 | ENSG00000233093 | 3.00 |
| RNU2-68P | ENSG00000222810 | 2.77 |

**Step 8 — PCA**

Data is scaled (zero mean, unit variance per gene) and then PCA is run on the 500 HVGs, computing 50 principal components. An elbow/variance plot is used to decide how many PCs capture most of the variance.

```
PCA shape: (252, 50)
```

**Step 9 — Neighborhood Graph + UMAP**

A k-nearest neighbor graph is built on the first 20 PCA components (k=15 neighbors). This graph is what both clustering and UMAP use.

UMAP then projects the neighborhood graph into 2D for visualization:

```
UMAP shape: (252, 2)
```

**Step 10 — Leiden Clustering**

Leiden is a community detection algorithm that finds groups of cells that are more connected to each other than to the rest of the dataset. Three resolutions were tested:

| Resolution | Clusters found |
|------------|---------------|
| 0.10 | 2 |
| 0.30 | 2 |
| 0.50 | 3 |

Resolution 0.3 was selected as the primary clustering since it gave a stable, biologically interpretable result with 2 clusters.

```
Cluster sizes (leiden, res=0.3):
  Cluster 0: 173 cells
  Cluster 1:  79 cells
```

**Steps 11–12 — QC re-assessment + Marker Genes**

After clustering, QC metrics are re-examined on the UMAP to check whether any clusters consist of low-quality cells. Marker genes per cluster are identified using the **Wilcoxon rank-sum test** (comparing each cluster against all others).

Top marker genes:

| Cluster | Top markers |
|---------|-------------|
| Cluster 0 | RPL3P12, RP11-558O12.1, RPL23AP83, RPSAP15, RPL19P21 |
| Cluster 1 | FTLP2, FTH1P8, AP1S2, CYBB, RAC1P4 |

**Steps 13–14 — Known marker genes + Cell type annotation**

Since the dataset is chrX-only, standard immune cell markers (CD3, CD19, etc.) aren't available — those genes are on autosomes. Instead, marker genes expressed on chrX were used. Available immune-relevant chrX genes found in this dataset include: `IL2RG`, `CXCR3`, `PIGA`, `TMSB4X`, `ELK1`, `RPS4X`, `RPL10`, `BCAP31`, `SEPT6`, `TAZ`, `FLNA`, `MECP2`.

Based on marker gene patterns:

| Cluster | Cell type | Cell count |
|---------|-----------|-----------|
| 0 | T/NK Lymphocytes | 173 |
| 1 | B Lymphocytes | 79 |

**Step 16 — Save results**

```
Output: pbmc_chrX_analyzed.h5ad (0.93 MB)

Contains:
  - Raw counts          → layers["counts"]
  - Log-normalized      → X and layers["log_normalized"]
  - PCA embedding       → obsm["X_pca"] (252 × 50)
  - UMAP embedding      → obsm["X_umap"] (252 × 2)
  - Leiden clusters     → obs["leiden"]
  - Cell type labels    → obs["cell_type"]
  - QC metrics          → obs columns
```

### Output files (in `results/`)

| File | Description |
|------|-------------|
| `pbmc_chrX_analyzed.h5ad` | Complete AnnData object with all embeddings, clusters, and annotations |
| `pbmc_cell_annotations.csv` | Flat CSV of all 252 cells × 19 metadata columns (QC metrics, doublet scores, cluster IDs, cell types) |

---

## Part 3 — AnnData Deep Dive

### What this part does

AnnData is the data structure that holds everything together throughout this pipeline — the count matrix, cell metadata, gene metadata, embeddings, clusters, and more. This notebook is a structured tutorial on how AnnData actually works: what each slot stores, how to access and subset it, and how to add new information.

**Input:** `pbmc_chrX_analyzed.h5ad` from Part 2

### The AnnData structure

When you load the Part 2 output:

```
AnnData object with n_obs × n_vars = 252 × 751
Shape: 252 cells × 751 genes
```

Here's what each slot contains in our specific dataset:

| Slot | What it stores | Our data |
|------|---------------|----------|
| `adata.X` | Active data matrix | Log-normalized counts — sparse CSR matrix, float32, 79.5% zeros |
| `adata.obs` | Cell-level metadata (DataFrame) | 19 columns: QC metrics, doublet scores, 4 leiden results, cell type |
| `adata.var` | Gene-level metadata (DataFrame) | 14 columns: Ensembl IDs, names, HVG flags, per-gene stats |
| `adata.obsm` | Multi-dimensional cell arrays | `X_pca` (252×50), `X_umap` (252×2) |
| `adata.layers` | Alternative data matrices | `counts` (raw int64), `log_normalized` (float32) |
| `adata.uns` | Unstructured metadata dict | 16 keys: colors, clustering params, PCA variance, neighbors config |
| `adata.obsp` | Cell×cell pair matrices | `connectivities` (252×252 float32), `distances` (252×252 float64) |

### Slot-by-slot exploration

**`adata.X` — The count matrix**

The active matrix holds log-normalized counts. Stored as a sparse matrix because most values are zero — a single cell only expresses a fraction of all genes.

```
Type:   scipy.sparse.csr_matrix
Shape:  (252, 751)
Non-zero entries:  38,710 out of 189,252 total
Sparsity:          79.5% zeros — only 20.5% of entries are non-zero
```

Example — first 5 cells × 5 genes (log-normalized values):
```
                  PLCXD1  GTPBP6  LINC00685  PPP2R3B  SHOX
AAAGAACCAATGGCAG    0.0   0.000      0.000      0.0   0.0
AAAGGATAGTAGACAT    0.0   0.528      0.000      0.0   0.0
AAAGGATCACCGGCTA    0.0   0.000      0.000      0.0   0.0
AAAGGATTCCGTTTCG    0.0   0.486      0.486      0.0   0.0
AAAGTCCCACCAGCCA    0.0   0.387      0.000      0.0   0.0
```

**`adata.obs` — Cell annotations**

A Pandas DataFrame indexed by the 16-character 10X cell barcode. Each row is one cell.

| Column | Type | Description |
|--------|------|-------------|
| `n_genes_by_counts` | int | Number of genes detected in this cell |
| `total_counts` | int | Total UMI count for this cell |
| `pct_counts_ribo` | float | % of counts from ribosomal genes |
| `doublet_score` | float | Scrublet score (0=clean, 1=likely doublet) |
| `predicted_doublet` | bool | Scrublet's binary doublet call |
| `leiden_res_0.10/0.30/0.50` | category | Cluster ID at each tested resolution |
| `leiden` | category | Final cluster assignment (res=0.3) |
| `cell_type` | str | Annotated label (T/NK Lymphocytes or B Lymphocytes) |
| `high_ribo` | bool | True if ribosomal content > 75th percentile |

Cell type distribution:
```
T/NK Lymphocytes    173 cells
B Lymphocytes        79 cells
```

Ribosomal content: 75th percentile threshold = 55.30% → 63 cells flagged as high-ribo, 189 cells normal.

**`adata.var` — Gene annotations**

A Pandas DataFrame indexed by gene name.

```
Shape: (751, 14)

Key columns:
  gene_ids          Ensembl ID (e.g., ENSG00000182378)
  gene_names        Gene symbol (e.g., PLCXD1)
  ribo              True if ribosomal gene
  n_cells_by_counts Number of cells where this gene is detected
  highly_variable   True for the 500 selected HVGs
  dispersions_norm  Normalized dispersion (HVG selection metric)
```

Example — metadata for gene `PLCXD1`:
```
gene_ids               ENSG00000182378
n_cells_by_counts      20
pct_dropout_by_counts  92.1%     ← only detected in 8% of cells
highly_variable        True
```

**`adata.obsm` — Embeddings**

Stores per-cell multi-dimensional arrays. All have 252 rows (one per cell).

```
X_pca:  (252, 50)  — 50 principal components
X_umap: (252, 2)   — 2D UMAP coordinates for visualization
```

PCA coordinates for first 5 cells (PCs 1–3):
```
                     PC1     PC2     PC3
AAAGAACCAATGGCAG -0.9848  3.0696  2.3816
AAAGGATAGTAGACAT -2.0650 -1.5160  0.5471
AAAGGATCACCGGCTA -1.9880 -1.5364 -0.8044
AAAGGATTCCGTTTCG  3.6262 -0.4807  0.9270
AAAGTCCCACCAGCCA -1.0071 -0.2777 -0.3916
```

UMAP coordinates for first 5 cells:
```
                      UMAP1    UMAP2
AAAGAACCAATGGCAG  19.3365   1.9986
AAAGGATAGTAGACAT  21.0020   6.2876
AAAGGATCACCGGCTA  19.8560   7.0598
AAAGGATTCCGTTTCG  -0.2017  12.0330
AAAGTCCCACCAGCCA  18.4466   5.9026
```

Cells in the same cluster have similar UMAP coordinates — that's the whole point.

**`adata.layers` — Multiple versions of the data**

Layers let you keep multiple versions of the count matrix simultaneously without duplicating the AnnData object. After Part 3:

```
counts:             (252, 751) int64   — raw integer UMI counts from Part 1
log_normalized:     (252, 751) float32 — log1p-normalized counts (same as X)
counts_per_million: (252, 751) float32 — CPM normalization (added in Part 3)
```

The transformation chain:
```
raw count (int) → library-size normalize → log1p → log-normalized value (float)
```

**`adata.uns` — Unstructured metadata**

A Python dictionary with 16 keys:

```
Keys: cell_type_colors, hvg, leiden_colors, leiden_res_0.10, leiden_res_0.10_colors,
      leiden_res_0.30, leiden_res_0.30_colors, leiden_res_0.50, leiden_res_0.50_colors,
      log1p, neighbors, pca, predicted_doublet_colors, rank_genes_groups, scrublet, umap
```

Project metadata added in Part 3:
```python
{
  'project_name': 'scRNA-seq 10X Preprocessing Project',
  'dataset':      '1k PBMC Healthy Donor (10X Genomics v3)',
  'genome':       'Human hg19 chrX',
  'part1_tool':   'RNA STARsolo + DropletUtils (Galaxy)',
  'part2_tool':   'Scanpy (Python)',
  'cells':        252,
  'genes':        751
}
```

**`adata.obsp` — Cell-cell pair matrices**

Stores the k-nearest neighbor graph computed during the neighbors step:

```
connectivities: (252, 252) float32 — weighted adjacency matrix
distances:      (252, 252) float64 — pairwise distances in PCA space
```

These sparse matrices encode which cells are neighbors of which. Cells in the same cluster have high connectivity values; cells from different clusters have near-zero values. Visualized as a heatmap, the two-cluster structure is clearly visible as two blocks of high connectivity along the diagonal.

### Key concepts demonstrated in Part 3

**Views vs Copies** — Subsetting AnnData (e.g., `adata[:10, :5]`) returns a *view* by default, which references the original data without copying it. Converting to a copy (`adata_subset.copy()`) creates an independent object. If you try to modify a view in-place, it automatically converts to a copy.

```
View:  is_view = True  → no extra memory, linked to parent
Copy:  is_view = False → independent, safe to modify
```

**Flexible subsetting** — AnnData supports slicing by:
- Integer index: `adata[:10, :5]`
- Name (barcode or gene): `adata[['AAAGAACCAATGGCAG', ...], :]`
- Boolean mask: `adata[adata.obs['cell_type'] == 'B Lymphocytes', :]` → 79 cells

**CPM normalization** — A new `counts_per_million` layer is added, showing how you can store multiple normalizations without duplicating the whole object.

**Export to DataFrames** — `adata.to_df()` converts `adata.X` to a Pandas DataFrame for inspection or CSV export. `adata.obs` and `adata.var` are already DataFrames and can be exported directly.

### Output file

```
pbmc_chrX_final.h5ad (1.07 MB)

Contains everything from Part 2, plus:
  - CPM layer (counts_per_million)
  - high_ribo flag in obs
  - Project metadata in uns
  - Verified read-back after writing
```

---

## Complete Pipeline Summary

```
Raw FASTQ reads (4 files, 2 lanes)
         │
         ▼  RNA STARsolo (Galaxy)
Aligned + demultiplexed: 5,200 barcodes
         │
         ▼  DropletUtils EmptyDrops (Galaxy)
252 high-quality cells × 2,392 chrX genes
  → barcodes.tsv / genes.tsv / matrix.mtx
         │
         ▼  Load into AnnData (Scanpy, Python)
252 × 2,392
         │
         ▼  QC + Gene filtering (min 3 cells)
252 × 751  (removed 1,641 low-expression genes)
         │
         ▼  Scrublet doublet detection
Doublet scores added per cell
         │
         ▼  Normalize (library-size) + log1p + HVG selection
500 highly variable genes identified
         │
         ▼  PCA (50 PCs) → k-NN graph (k=15) → UMAP
2D visualization ready
         │
         ▼  Leiden clustering (res=0.3)
2 clusters: 173 cells | 79 cells
         │
         ▼  Wilcoxon marker genes + manual annotation
T/NK Lymphocytes (173) | B Lymphocytes (79)
         │
         ▼  AnnData enrichment (Part 3)
CPM layer + project metadata + obsp exploration
         │
         ▼
pbmc_chrX_final.h5ad (1.07 MB)
```

---

## How to Run

### Part 1 (Galaxy — no installation required)
1. Create a free account at [usegalaxy.org](https://usegalaxy.org)
2. Upload the 6 input files from [Zenodo 3457880](https://zenodo.org/record/3457880)
3. Follow the [Galaxy GTN tutorial](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html)
4. Run STARsolo → MultiQC → DropletUtils and download the 3 output files

### Parts 2 & 3 (Python / Google Colab)
1. Open the notebooks in Google Colab or Jupyter
2. Run Step 0 once to install dependencies:
   ```bash
   pip install scanpy anndata scrublet
   ```
3. Upload the 3 MTX files from Part 1 when prompted in Part 2, or the `.h5ad` file for Part 3
4. Run all cells sequentially

---

## Tools & Versions

| Tool | Version | Part | Platform |
|------|---------|------|----------|
| RNA STARsolo | 2.7.11a | 1 | Galaxy |
| DropletUtils | 1.10.0 | 1 | Galaxy |
| MultiQC | 1.27 | 1 | Galaxy |
| Scanpy | 1.12.1 | 2, 3 | Python |
| AnnData | 0.12.10 | 2, 3 | Python |
| NumPy | 2.0.2 | 2, 3 | Python |
| Pandas | 2.3.3 | 2, 3 | Python |
| Scrublet | latest | 2 | Python |

---

## References

- Wolf FA, Angerer P, Theis FJ (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 19:15.
- Lun ATL et al. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. *Genome Biology* 20:63.
- Wolock SL, Lopez R, Klein AM (2019). Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Systems* 8(4):281–291.
- Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports* 9:5233.
- 10X Genomics dataset: [1k PBMCs from a Healthy Donor (v3 chemistry)](https://www.10xgenomics.com/resources/datasets)
- Galaxy GTN tutorial: [Pre-processing of 10X Single-Cell RNA Datasets](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html)
- scverse clustering tutorial: [https://scverse-tutorials.readthedocs.io](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/clustering.html)
- AnnData documentation: [https://anndata.readthedocs.io](https://anndata.readthedocs.io)

---

