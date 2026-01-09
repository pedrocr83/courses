# Environment Setup

## R Environment (Primary for CellChat, NicheNet)

### Install Core Packages

```r
# Install devtools if needed
install.packages("devtools")

# CellChat
devtools::install_github("sqjin/CellChat")

# NicheNet
devtools::install_github("saeyslab/nichenetr")

# LIANA (R version)
devtools::install_github("saezlab/liana")

# Bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "SingleCellExperiment",
    "ComplexHeatmap",
    "Seurat"
))

# Visualization
install.packages(c(
    "ggplot2",
    "igraph",
    "circlize",
    "ggalluvial",
    "patchwork"
))
```

### Verify R Installation

```r
library(CellChat)
library(nichenetr)
library(liana)
library(Seurat)

cat("CellChat version:", packageVersion("CellChat"), "\n")
cat("NicheNet version:", packageVersion("nichenetr"), "\n")
cat("LIANA version:", packageVersion("liana"), "\n")
```

---

## Python Environment (CellPhoneDB, LIANA-py)

### Conda Setup

```bash
# Create environment
conda create -n ccc-course python=3.11 -y
conda activate ccc-course

# Core packages
pip install scanpy anndata pandas numpy matplotlib seaborn

# CellPhoneDB
pip install cellphonedb

# LIANA Python
pip install liana

# Jupyter
pip install jupyterlab
```

### Verify Python Installation

```python
import scanpy as sc
import cellphonedb
import liana

print(f"scanpy: {sc.__version__}")
print("CellPhoneDB installed")
print("LIANA installed")
```

---

## Demo Data

### PBMC 3k (Primary dataset)

#### Python
```python
import scanpy as sc
adata = sc.datasets.pbmc3k_processed()
```

#### R (Seurat)
```r
library(Seurat)
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k")
```

### CellChat Tutorial Data

```r
# Load CellChat example data
library(CellChat)
data(Skin_Merged)  # Human skin dataset
```

---

## Ligand-Receptor Databases

### CellChatDB

```r
library(CellChat)
CellChatDB <- CellChatDB.human  # or CellChatDB.mouse
showDatabaseCategory(CellChatDB)
```

### CellPhoneDB

Automatically downloaded when running analysis.

### OmniPath (via LIANA)

```r
library(liana)
# LIANA automatically accesses OmniPath
```

---

## Directory Structure

```
cell_cell_communication_course/
├── data/
│   ├── raw/
│   └── processed/
├── labs/
├── assignments/
├── results/
│   ├── cellchat/
│   ├── cellphonedb/
│   ├── nichenet/
│   └── liana/
├── figures/
├── scripts/
└── resources/
```

Create this structure:

```bash
mkdir -p data/{raw,processed} results/{cellchat,cellphonedb,nichenet,liana} figures scripts
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| CellChat install fails | Install dependencies first: `BiocManager::install("ComplexHeatmap")` |
| NicheNet download slow | Pre-download networks from GitHub |
| CellPhoneDB memory error | Subsample cells or use statistical analysis only |
| LIANA version conflict | Create fresh environment |
| Seurat v5 compatibility | Check CellChat version for Seurat v5 support |

---

## Hardware Requirements

| Analysis | RAM | Time (3k cells) |
|----------|-----|-----------------|
| CellPhoneDB | 8 GB | 10-30 min |
| CellChat | 8 GB | 5-15 min |
| NicheNet | 16 GB | 15-45 min |
| LIANA (all methods) | 16 GB | 20-60 min |

---

## Recommended Resources

### Papers to Download
- Armingol et al., 2021 (CCC review)
- Jin et al., 2021 (CellChat)
- Browaeys et al., 2020 (NicheNet)
- Efremova et al., 2020 (CellPhoneDB)

### Tutorials
- CellChat: https://sqjin.github.io/CellChat/
- NicheNet: https://github.com/saeyslab/nichenetr
- CellPhoneDB: https://www.cellphonedb.org/
- LIANA: https://saezlab.github.io/liana/

