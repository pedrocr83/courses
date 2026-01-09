# Environment Setup

## Python Environment (scvi-tools)

```bash
conda create -n integration-course python=3.11 -y
conda activate integration-course

# Core
pip install scanpy anndata pandas numpy matplotlib seaborn

# Integration methods
pip install harmonypy  # Harmony
pip install scvi-tools  # scVI, scANVI

# Benchmarking
pip install scib  # Single-cell integration benchmarking

# Jupyter
pip install jupyterlab
```

### GPU Support (Optional, for scVI)

```bash
# If you have CUDA
pip install scvi-tools[cuda]
```

### Verify Installation

```python
import scanpy as sc
import harmonypy
import scvi
import scib
print(f"scanpy: {sc.__version__}")
print(f"scvi: {scvi.__version__}")
```

---

## R Environment (Seurat, Harmony)

```r
# Core
install.packages("Seurat")

# Integration
install.packages("harmony")
BiocManager::install("batchelor")  # MNN

# Evaluation
# devtools::install_github("theislab/kBET")

# Visualization
install.packages(c("ggplot2", "patchwork", "cowplot"))
```

### Verify R Installation

```r
library(Seurat)
library(harmony)
library(batchelor)
packageVersion("Seurat")
```

---

## Demo Datasets

### Pancreas Integration Benchmark

Classic multi-technology benchmark dataset.

```python
# Download from scib tutorials
import scanpy as sc
# See scib documentation for download
```

### Multi-Donor PBMC

```python
# Multiple PBMC donors with batch effects
# Available from cellxgene
```

---

## Benchmarking Metrics (scIB)

```python
import scib

# Key metrics:
# - iLISI: Integration LISI (batch mixing)
# - cLISI: Cell-type LISI (biological conservation)
# - kBET: k-nearest neighbor batch effect test
# - ARI: Adjusted Rand Index
# - NMI: Normalized Mutual Information
```

---

## Directory Structure

```
data_integration_course/
├── data/
│   ├── batch1/
│   ├── batch2/
│   └── integrated/
├── labs/
├── assignments/
├── results/
│   ├── harmony/
│   ├── seurat/
│   ├── scvi/
│   └── benchmarks/
├── figures/
└── resources/
```

---

## Memory Requirements

| Method | RAM (10k cells) | GPU |
|--------|-----------------|-----|
| Harmony | 8 GB | No |
| Seurat CCA | 16 GB | No |
| MNN | 8 GB | No |
| scVI | 8 GB | Optional |
| scANVI | 8 GB | Optional |

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| scVI slow | Use GPU or reduce epochs |
| Memory error | Subsample or use sparse matrices |
| Harmony fails | Check input is PCA embedding |
| kBET slow | Use subset of cells |

