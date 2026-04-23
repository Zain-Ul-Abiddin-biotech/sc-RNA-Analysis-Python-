# 🧬 scRNA-seq T Cell Subtype Identification — Python / Scanpy

> **Who is this for?**  
> Python beginners who want to understand and run a complete single-cell RNA-seq
> analysis. Every section explains *what* the code does, *why* it matters, and
> *what you can change* for your own data.

---

## 📋 Table of Contents

1. [What does this notebook do?](#1-what-does-this-notebook-do)
2. [Requirements](#2-requirements--what-to-install)
3. [The Dataset](#3-the-dataset)
4. [Section-by-Section Walkthrough](#4-section-by-section-walkthrough)
5. [All Bugs Fixed](#5-all-bugs-fixed)
6. [Tweaking the Pipeline](#6-tweaking-the-pipeline)
7. [Common Errors & Fixes](#7-common-errors--fixes)
8. [Project Structure](#8-project-structure)
9. [Further Reading](#9-further-reading)

---

## 1. What does this notebook do?

This notebook runs a full **single-cell RNA sequencing (scRNA-seq)** analysis on
human CD4+ T cells, identifying four biologically distinct subtypes.

In plain English: you have ~5,000 individual cells. For each cell you know how
active each of ~14,700 genes is. This pipeline:

- **Cleans** the data (removes dead cells, doublets, empty droplets)
- **Normalises** expression so cells are comparable
- **Clusters** cells that look similar together
- **Visualises** clusters on a 2D UMAP plot
- **Identifies** what T cell subtype each cluster is

**Final output:** a UMAP where every dot is one cell, coloured by type:

| Cell type | Full name | Biology |
|-----------|-----------|---------|
| TN | T Naive | Resting; never seen antigen |
| TCM | T Central Memory | Long-lived; recirculate in blood |
| TEM | T Effector Memory | Tissue-resident; rapid responders |
| TEMRA | Terminally differentiated RA | Cytotoxic; end-stage |

---

## 2. Requirements — what to install

**Python 3.9 or higher required.**

```bash
pip install scanpy anndata numpy pandas matplotlib
```

Or with conda (recommended):

```bash
conda create -n scrna python=3.10
conda activate scrna
pip install scanpy anndata numpy pandas matplotlib
```

Verify the install:
```python
import scanpy as sc
print(sc.__version__)   # Should be 1.9.0 or higher
```

---

## 3. The Dataset

This pipeline uses `scdr.h5ad` — human CD4+ T cells sequenced with 10x Genomics.

### What is an `.h5ad` file?

The standard single-cell data format in Python. One file stores everything:

```
scdr.h5ad
├── .X       expression matrix (cells x genes)
├── .obs     cell metadata: cell.type, cytokine.condition, donor.id, batch, etc.
├── .var     gene metadata
└── .obsm    dimensionality reductions (added by the pipeline)
```

| Property | Value |
|----------|-------|
| Cells | ~5,000 CD4+ T cells |
| Genes | ~14,700 |
| Cell types | TN, TCM, TEM, TEMRA |

**To use your own data:**
```python
scdr = sc.read_h5ad('your_file.h5ad')

# 10x Genomics folder format:
scdr = sc.read_10x_mtx('path/to/filtered_feature_bc_matrix/',
                        var_names='gene_symbols')
```

Public datasets: [CELLxGENE](https://cellxgene.cziscience.com/) · [Human Cell Atlas](https://www.humancellatlas.org/) · [GEO](https://www.ncbi.nlm.nih.gov/geo/)

---

## 4. Section-by-Section Walkthrough

---

### Section 1 — Load Libraries

```python
import scanpy as sc
sc.settings.verbosity = 2
sc.settings.seed = 0
np.random.seed(0)
```

| Library | Purpose |
|---------|---------|
| `scanpy` | Core single-cell analysis toolkit |
| `anndata` | Reads/writes `.h5ad` files |
| `numpy` | Fast array maths |
| `pandas` | DataFrames and tables |
| `matplotlib` | Plotting (used internally by scanpy) |

**Why set a random seed?**  
Leiden clustering and UMAP both use random initialisation. Without a fixed seed,
cluster numbers change every run — making Section 10's cell-type assignments
wrong on re-runs. `seed=0` makes results reproducible.

---

### Section 2 — Load Data & Explore

```python
scdr = sc.read_h5ad('scdr.h5ad')
print(scdr)   # n_obs = cells, n_vars = genes
```

Useful commands to explore:
```python
scdr.obs.columns.tolist()             # All cell metadata columns
scdr.obs['cell.type'].value_counts()  # Cells per type
scdr.obs.head()                       # First 5 rows
```

---

### Section 3 — Subset to Naive T Cells

```python
ids  = scdr.obs['cell.type'] == 'Naive'
scdr = scdr[ids, :].copy()   # .copy() is required — avoids view warnings
```

**Index notation:**
```python
scdr[row_selector, column_selector]
scdr[ids, :]   # keep selected cells, keep ALL genes
scdr[:, ids]   # keep ALL cells, keep selected genes
```

To analyse a different cell type:
```python
ids = scdr.obs['cell.type'] == 'Memory'   # change here
```

---

### Section 4 — Quality Control

Three types of bad cells to remove:

| Problem | Signal | What it is |
|---------|--------|-----------|
| Empty droplet | Very low `n_genes_by_counts` | No cell was captured |
| Doublet | Very high `n_genes_by_counts` | Two cells merged into one |
| Dying cell | High `pct_counts_mt` | Cytoplasm leaked; only mitochondria remain |

**Mitochondrial genes** start with `MT-` in humans (use `mt-` for mouse).
High mitochondrial % = likely dead or damaged cell.

**QC thresholds to adjust based on your violin plots:**

| Threshold | Default | When to change |
|-----------|---------|----------------|
| `n_genes < 2000` | 2000 | Raise for complex cells (neurons) |
| `pct_mt < 5` | 5 % | Raise to 10-25 % for heart/muscle tissue |

---

### Section 5 — Preprocessing

**Always in this order:**
```
Raw counts
  -> normalize_total()   each cell sums to 10,000
  -> log1p()             log(x + 1) applied to every value
Normalised log counts
```

**Why normalise?** Some cells are sequenced deeper — they look artificially
"higher" for every gene. Rescaling to 10,000 makes cells comparable.

**Why log-transform?** A gene 10x more expressed = 9,000 count difference.
After log1p, same 10x change = ~2.3. Much easier to work with statistically.

---

### Section 6 — Highly Variable Genes & Scaling

**Why select HVGs?** Most genes are equally expressed in all cells and carry
no information. The most variable genes (different between cell types) are
what we care about. Keeping ~2,000 HVGs speeds up all downstream steps.

**The `.raw` snapshot:**
```python
scdr.raw = scdr.copy()   # save normalised log counts before scaling
```
Required later for differential expression and UMAP gene-expression colouring.

**Scaling:** centers each gene to mean=0, std=1. Without this, highly expressed
genes dominate PCA just because of their large numbers.

---

### Section 7 — PCA

PCA compresses 2,000 HVG dimensions into ranked principal components (PCs).

**How to read the scree plot:** find the "elbow" where the curve flattens.
Use that PC number as `n_pcs` in Section 8.

**Scanpy module guide:**
| Module | Job |
|--------|-----|
| `sc.pp` | preprocessing (filtering, normalising, scaling) |
| `sc.tl` | tools (PCA, clustering, UMAP, DE testing) |
| `sc.pl` | all visualisations |

---

### Section 8 — Clustering & UMAP

```python
sc.pp.neighbors(scdr, n_neighbors=10, n_pcs=30)  # build cell graph
sc.tl.leiden(scdr, resolution=1.0)                # find clusters
sc.tl.umap(scdr)                                  # project to 2D
```

| Parameter | Default | Effect |
|-----------|---------|--------|
| `n_pcs` | `30` | Use elbow plot value |
| `n_neighbors` | `10` | Higher = smoother UMAP, less local detail |
| `resolution` | `1.0` | Higher = more clusters; lower = fewer |

**UMAP caveats:**
- Cells close together = genuinely similar (reliable)
- Distance *between* clusters = not always meaningful

---

### Section 9 — Marker Genes

```python
sc.tl.rank_genes_groups(scdr, 'leiden', method='wilcoxon')
```

For each cluster, finds genes most specifically expressed there.
Look up the top 5-10 genes per cluster to identify the cell type.

**Known T cell markers:**

| Marker genes | Cell type |
|---|---|
| `CCR7`, `SELL`, `TCF7`, `LEF1` | T Naive (TN) |
| `CCR7`, `IL7R`, `S100A4` | T Central Memory (TCM) |
| `GZMK`, `CCL5`, `EOMES` | T Effector Memory (TEM) |
| `GZMB`, `PRF1`, `KLRG1`, `CX3CR1` | TEMRA |

Lookup databases: [GeneCards](https://www.genecards.org/) · [CellMarker 2.0](http://biocc.hrbmu.edu.cn/CellMarker/) · [PanglaoDB](https://panglaodb.se/)

---

### Section 10 — Cell Type Annotation

```python
scdr.obs.loc[scdr.obs.leiden.isin(['0','3','5']), 'cell_identity'] = 'TN'
scdr.obs.loc[scdr.obs.leiden.isin(['1','2']),     'cell_identity'] = 'TCM'
scdr.obs.loc[scdr.obs.leiden == '4',              'cell_identity'] = 'TEM'
scdr.obs.loc[scdr.obs.leiden == '6',              'cell_identity'] = 'TEMRA'
```

> **Warning:** Cluster numbers are not stable. They change if you alter
> `resolution`, `n_neighbors`, or `n_pcs`. Always re-check marker plots
> before re-using these assignments.

---

### Section 11 — Save Results

```python
scdr.write_h5ad('scdr_cell_classified.h5ad')

# Reload later:
scdr = sc.read_h5ad('scdr_cell_classified.h5ad')
```

---

## 5. All Bugs Fixed

| # | Original cell | Bug | Fix applied |
|---|--------------|-----|-------------|
| 1 | Cell 2 | `np.flip(...)[range(10)]` — list index deprecated in NumPy >= 1.24 | Changed to `[:10]` slice |
| 2 | Cells 4-7 | Filtered to Naive cells (cells 4-6) then **reloaded full dataset** (cell 7), silently discarding all work | Consolidated: reload once, subset, proceed |
| 3 | Cell 6 | `scdr.var = pd.DataFrame(...)` replaced the **entire** var table, wiping gene names needed by QC | Changed to `scdr.var['num_reads'] = ...` (add column only) |
| 4 | Cells 17, 19 | `scdr.X.todense()` deprecated in SciPy >= 1.14 / NumPy >= 2.0 | Replaced with `.toarray()` plus a dense fallback |
| 5 | Cell 31 | `sc.tl.leiden()` without `resolution` — deprecated in Scanpy >= 1.9 | Added `resolution=1.0` explicitly |
| 6 | Cell 33 | `method='t-test'` assumes normality — invalid for scRNA-seq data | Changed to `method='wilcoxon'` (community standard) |
| 7 | Cells 34-35 | Chained `|` conditions for cluster assignment — verbose and error-prone | Replaced with `.isin([...])` |
| 8 | Everywhere | No random seed — cluster numbers change between runs | Added `sc.settings.seed = 0` and `np.random.seed(0)` |
| 9 | Multiple cells | Stray `\\n"` artifacts from Jupyter JSON export in comment strings | Cleaned up throughout |
| 10 | Cell 36 | Bare `scdr` with no useful output at the end | Replaced with `print(scdr)` for an explicit final summary |

---

## 6. Tweaking the Pipeline

| What to change | Parameter | Default |
|----------------|-----------|---------|
| Input file | `sc.read_h5ad('...')` | `'scdr.h5ad'` |
| Cell type to analyse | `== 'Naive'` | `'Naive'` |
| Min genes per cell | `min_genes=200` | `200` |
| Doublet gene cutoff | `< 2000` | `2000` |
| Mitochondrial cutoff | `< 5` | `5 %` |
| Mouse mito prefix | `'MT-'` | Change to `'mt-'` |
| Number of PCs | `n_pcs=30` | `30` (use scree plot) |
| Neighbour count | `n_neighbors=10` | `10` |
| Cluster resolution | `resolution=1.0` | `1.0` |
| DE test method | `method='wilcoxon'` | `'wilcoxon'` |
| Cluster-to-celltype map | `.isin([...])` blocks | Based on your marker genes |

---

## 7. Common Errors & Fixes

| Error | Cause | Fix |
|-------|-------|-----|
| `ModuleNotFoundError: scanpy` | Not installed | `pip install scanpy` |
| `FileNotFoundError: scdr.h5ad` | Wrong path | Use `os.path.abspath('scdr.h5ad')` to find full path |
| `ImplicitModificationWarning` | Modified a view | Add `.copy()` after every `scdr[...]` slice |
| UMAP looks like one blob | Too few PCs or over-filtering | Increase `n_pcs`; loosen QC thresholds |
| All cells in one cluster | Resolution too low | Increase `resolution` in `sc.tl.leiden()` |
| 50+ tiny clusters | Resolution too high | Decrease `resolution` |
| `UserWarning: log1p base` | Known Scanpy reload bug | Already fixed: `scdr.uns['log1p'] = {'base': None}` |
| Unassigned cells after annotation | Leiden produced new cluster numbers | Re-run Section 9 and re-map clusters |

---

## 8. Project Structure

```
your-project/
│
├── scRNA_SEQ_ANALYSIS_corrected.ipynb   <- This notebook
├── README.md                            <- This file
│
└── data/
    ├── scdr.h5ad                        <- Raw input (download separately)
    └── scdr_preprocessed.h5ad          <- Generated by the notebook
```

> **Do not commit `.h5ad` files to GitHub** — they are typically 100 MB+.  
> Add to `.gitignore`:
> ```
> *.h5ad
> data/
> ```

---

## 9. Further Reading

- [Scanpy tutorials](https://scanpy.readthedocs.io/en/stable/tutorials.html) — official step-by-step guides
- [AnnData documentation](https://anndata.readthedocs.io/) — understanding the data structure
- [CELLxGENE](https://cellxgene.cziscience.com/) — browse and download public datasets
- [CellMarker 2.0](http://biocc.hrbmu.edu.cn/CellMarker/) — cell-type marker lookup
- [Luecken & Theis (2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746) — best-practices paper
- [Scanpy vs Seurat benchmark](https://www.nature.com/articles/s41592-019-0692-4) — Nature Methods

---

*Corrected from original Jupyter notebook — all 10 bugs documented in Section 5.*
