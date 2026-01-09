# Environment Setup

## Python Environment (Scanpy - Recommended)

```bash
conda create -n clustering-course python=3.11 -y
conda activate clustering-course

# Core packages
pip install scanpy anndata pandas numpy matplotlib seaborn

# Clustering
pip install leidenalg python-igraph

# Automated annotation
pip install celltypist

# Jupyter
pip install jupyterlab
```

### Verify Python Installation

```python
import scanpy as sc
import celltypist
print(f"scanpy: {sc.__version__}")
print("celltypist installed")
```

---

## R Environment (Seurat)

```r
# Core
install.packages("Seurat")
install.packages(c("ggplot2", "dplyr", "patchwork"))

# Annotation
BiocManager::install(c("SingleR", "celldex"))
BiocManager::install("scmap")

# Visualization
install.packages(c("pheatmap", "clustree"))
```

### Verify R Installation

```r
library(Seurat)
library(SingleR)
packageVersion("Seurat")
```

---

## Demo Datasets

### PBMC 3k (Primary)

```python
import scanpy as sc
adata = sc.datasets.pbmc3k()
```

```r
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k")
```

### PBMC 10k (Assignments)

Download from 10x Genomics website.

---

## Marker Gene Resources

### PanglaoDB
- Website: https://panglaodb.se/
- Download marker lists for cell types

### CellMarker 2.0
- Website: http://bio-bigdata.hrbmu.edu.cn/CellMarker/

### CellTypist Models
```python
import celltypist
celltypist.models.download_models()
celltypist.models.models_description()
```

---

## Directory Structure

```
clustering_annotation_course/
├── data/
│   ├── raw/
│   └── processed/
├── labs/
├── assignments/
├── results/
│   ├── clustering/
│   └── annotation/
├── figures/
└── resources/
    └── markers/
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Leiden fails | `pip install leidenalg python-igraph` |
| SingleR slow | Use subset of reference |
| Memory errors | Subsample cells for testing |

