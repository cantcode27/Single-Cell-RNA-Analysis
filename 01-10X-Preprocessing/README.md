# Part 1: Pre-processing of 10X Single-Cell RNA Datasets

## Overview
This section contains the Galaxy-based preprocessing pipeline that converts raw 10X Genomics FASTQ files into a clean, filtered count matrix using RNA STARsolo and DropletUtils.

## Tools Used
| Tool | Version | Platform |
|------|---------|----------|
| RNA STARsolo | 2.7.11a | Galaxy (usegalaxy.org) |
| DropletUtils | 1.10.0 | Galaxy |
| MultiQC | 1.27 | Galaxy |

## Input Data
Downloaded from [Zenodo record 3457880](https://zenodo.org/record/3457880):
- `subset_pbmc_1k_v3_S1_L001_R1_001.fastq.gz` — Lane 1 barcodes
- `subset_pbmc_1k_v3_S1_L001_R2_001.fastq.gz` — Lane 1 cDNA
- `subset_pbmc_1k_v3_S1_L002_R1_001.fastq.gz` — Lane 2 barcodes
- `subset_pbmc_1k_v3_S1_L002_R2_001.fastq.gz` — Lane 2 cDNA
- `3M-february-2018.txt.gz` — 10X v3 barcode whitelist
- `Homo_sapiens.GRCh37.75.gtf` — Gene annotation

## Pipeline Steps

### Step 1: RNA STARsolo
- Aligns reads to hg19 chrX reference
- Demultiplexes cells using barcode whitelist
- Quantifies gene expression per cell

### Step 2: MultiQC
- Assesses mapping quality from STARsolo log

### Step 3: DropletUtils — EmptyDrops
- Filters real cells from empty droplets
- Parameters: lower=200, FDR=0.01
- **Result: 252 high-quality cells**

## Output Files (in `results/`)
| File | Description |
|------|-------------|
| `barcodes.tsv` | 252 cell barcodes |
| `genes.tsv` | 2,392 gene IDs and names (chrX) |
| `matrix.mtx` | Sparse count matrix (MatrixMarket format) |

## Key Metrics
- **Cells detected by STARsolo**: 5,200
- **High-quality cells (EmptyDrops)**: 252
- **Genes quantified**: 2,392 (chrX only)
- **Non-zero entries in matrix**: 39,002

## How to Reproduce
1. Go to [usegalaxy.org](https://usegalaxy.org)
2. Create a new history named `scRNA-seq 10X Preprocessing Project`
3. Upload the 6 input files from Zenodo
4. Follow the [Galaxy GTN tutorial](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html)
5. Run STARsolo → MultiQC → DropletUtils as described in METHODOLOGY.md

## Reference
- Galaxy Training Network: [Pre-processing of 10X Single-Cell RNA Datasets](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html)
