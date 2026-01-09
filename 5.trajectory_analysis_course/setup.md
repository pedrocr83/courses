# Environment Setup

## Python Environment (scVelo, CellRank)

```bash
conda create -n trajectory-course python=3.11 -y
conda activate trajectory-course

# Core
pip install scanpy anndata pandas numpy matplotlib seaborn

# Trajectory analysis
pip install scvelo  # RNA velocity
pip install cellrank  # Fate probabilities

# Loom file support (for velocity)
pip install loompy

# Jupyter
pip install jupyterlab
```

### Verify Python Installation

```python
import scanpy as sc
import scvelo as scv
import cellrank as cr

print(f"scanpy: {sc.__version__}")
print(f"scvelo: {scv.__version__}")
print(f"cellrank: {cr.__version__}")
```

---

## R Environment (Monocle3, Slingshot)

```r
# Monocle3
devtools::install_github("cole-trapnell-lab/monocle3")

# Slingshot
BiocManager::install("slingshot")

# Trajectory DE
BiocManager::install("tradeSeq")

# Visualization
install.packages(c("ggplot2", "viridis", "pheatmap"))
```

### Verify R Installation

```r
library(monocle3)
library(slingshot)
library(tradeSeq)
packageVersion("monocle3")
```

---

## Velocity Data Preparation

RNA velocity requires spliced/unspliced counts. Two options:

### Option 1: velocyto (Command line)

```bash
pip install velocyto

# Run on BAM file
velocyto run -b filtered_barcodes.tsv -o output/ \
    possorted_genome_bam.bam annotation.gtf
```

### Option 2: kallisto | bustools (Faster)

```bash
# See kb-python documentation
pip install kb-python

kb count --workflow lamanno \
    -i index.idx -g t2g.txt \
    -x 10xv3 -o output/ \
    R1.fastq.gz R2.fastq.gz
```

---

## Demo Datasets

### Pancreas Endocrine Differentiation

Classic trajectory dataset.

```python
import scvelo as scv
adata = scv.datasets.pancreas()
```

### Dentate Gyrus (with velocity)

```python
adata = scv.datasets.dentategyrus()
```

### Bone Marrow (Hematopoiesis)

```python
# Multiple lineages
adata = scv.datasets.bonemarrow()
```

---

## Directory Structure

```
trajectory_analysis_course/
├── data/
│   ├── raw/
│   ├── velocity/
│   └── processed/
├── labs/
├── assignments/
├── results/
│   ├── pseudotime/
│   ├── velocity/
│   └── cellrank/
├── figures/
└── resources/
```

---

## Computational Requirements

| Analysis | RAM | Time (5k cells) |
|----------|-----|-----------------|
| Diffusion PT | 8 GB | 2 min |
| Monocle3 | 8 GB | 5 min |
| Slingshot | 4 GB | 2 min |
| scVelo stochastic | 8 GB | 5 min |
| scVelo dynamical | 16 GB | 30 min |
| CellRank | 16 GB | 10 min |

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| No velocity data | Need to process with velocyto/kb |
| Monocle3 install fails | Check R version compatibility |
| CellRank memory | Use sparse mode |
| Velocity arrows chaotic | Try dynamical model |
| No clear trajectory | Data may not have continuous process |

