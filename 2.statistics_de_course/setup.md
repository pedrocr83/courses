# Environment Setup

## R Environment (Primary - Recommended)

### Option 1: Install Packages Directly

```r
# Install Bioconductor manager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Core DE packages
BiocManager::install(c(
    "DESeq2",
    "edgeR",
    "limma"
))

# Example datasets
BiocManager::install(c(
    "airway",
    "pasilla"
))

# Annotation and enrichment
BiocManager::install(c(
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "clusterProfiler",
    "enrichplot",
    "DOSE"
))

# Visualization
install.packages(c(
    "ggplot2",
    "pheatmap",
    "RColorBrewer",
    "ggrepel"
))

# Install EnhancedVolcano
BiocManager::install("EnhancedVolcano")

# Tidyverse for data wrangling
install.packages("tidyverse")
```

### Option 2: Conda Environment

```bash
# Create environment
conda create -n de-course r-base=4.3 -y
conda activate de-course

# Install Bioconductor packages via conda
conda install -c bioconda bioconductor-deseq2 bioconductor-edger bioconductor-limma
conda install -c bioconda bioconductor-airway bioconductor-clusterprofiler
conda install -c conda-forge r-ggplot2 r-pheatmap r-tidyverse
```

---

## Python Environment (Alternative)

```bash
# Create environment
conda create -n de-course-py python=3.11 -y
conda activate de-course-py

# Core packages
pip install scanpy anndata pandas numpy matplotlib seaborn

# DE analysis
pip install pydeseq2 diffxpy

# Enrichment
pip install gseapy

# Jupyter
pip install jupyterlab
```

---

## RStudio Setup

1. Download RStudio: https://posit.co/download/rstudio-desktop/
2. Set working directory to course folder
3. Create R Project for organization

---

## Verify Installation

### R

```r
# Check DESeq2
library(DESeq2)
packageVersion("DESeq2")

# Check edgeR
library(edgeR)
packageVersion("edgeR")

# Check airway dataset
library(airway)
data(airway)
airway
```

### Python

```python
import scanpy as sc
import pydeseq2
print(f"scanpy: {sc.__version__}")
```

---

## Demo Datasets

### Airway Dataset (Bulk RNA-seq)

Primary dataset for bulk DE labs. Already included in Bioconductor.

```r
library(airway)
data(airway)

# Explore
colData(airway)
head(assay(airway))
```

**Study:** Airway smooth muscle cells treated with dexamethasone
- 4 cell lines
- 2 conditions: treated vs untreated
- Paired design

### Pasilla Dataset (Bulk RNA-seq)

Alternative dataset for assignments.

```r
BiocManager::install("pasilla")
library(pasilla)
```

### PBMC 3k (scRNA-seq)

For single-cell DE labs.

```python
import scanpy as sc
adata = sc.datasets.pbmc3k()
```

---

## Directory Structure

```
statistics_de_course/
├── labs/
├── assignments/
├── data/
│   ├── raw/
│   └── processed/
├── results/
│   ├── deseq2/
│   ├── edger/
│   └── figures/
├── scripts/
└── resources/
```

Create this structure:

```bash
mkdir -p data/{raw,processed} results/{deseq2,edger,figures} scripts
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| DESeq2 install fails | Update Bioconductor: `BiocManager::install(version = "3.18")` |
| Memory errors | Increase R memory or use smaller dataset |
| Package conflicts | Create fresh R environment |
| Missing annotation DB | Install species-specific: `org.Hs.eg.db` for human |

---

## Recommended Resources

### Books (Free Online)
- Modern Statistics for Modern Biology (Holmes & Huber)
- RNA-seq Data Analysis (Law et al.)

### Vignettes
- DESeq2: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- edgeR: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

