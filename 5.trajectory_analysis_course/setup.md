# Environment Setup

> **Updated February 2026** — CellRank 2 protocol environment is now the primary setup.

---

## Python Environment (Primary: CellRank Protocol)

### Option A: Recommended — CellRank Protocol Environment

This setup follows the official [CellRank Protocol](https://github.com/theislab/cellrank_protocol) (Weiler & Theis, 2026).

```bash
# Create environment
conda create -n cellrank-course python=3.11 --yes
conda activate cellrank-course

# Install CellRank (includes all dependencies)
conda install -c conda-forge cellrank

# Clone the protocol repo (reference notebooks + data utilities)
git clone https://github.com/theislab/cellrank_protocol.git
cd cellrank_protocol
pip install -e ".[jupyter]"
cd ..

# Additional tools for the course
pip install scvelo       # RNA velocity
pip install loompy       # Velocity data format
pip install jupyterlab   # Notebook interface

# Register Jupyter kernel
python -m ipykernel install --user --name cellrank-course --display-name "cellrank-course"
```

### Option B: Minimal — pip only

```bash
conda create -n cellrank-course python=3.11 -y
conda activate cellrank-course

# Core
pip install scanpy anndata pandas numpy matplotlib seaborn

# Trajectory analysis
pip install scvelo       # RNA velocity
pip install cellrank     # Fate mapping (CellRank 2)
pip install loompy       # Velocity data format

# Jupyter
pip install jupyterlab
```

### Verify Python Installation

```python
import scanpy as sc
import scvelo as scv
import cellrank as cr

print(f"scanpy:   {sc.__version__}")
print(f"scvelo:   {scv.__version__}")
print(f"cellrank: {cr.__version__}")

# Verify CellRank 2 kernels are available
from cellrank.kernels import VelocityKernel, CytoTRACEKernel, PseudotimeKernel, RealTimeKernel
from cellrank.estimators import GPCCA
print("\nAll CellRank 2 kernels available!")
print("  - VelocityKernel")
print("  - CytoTRACEKernel")
print("  - PseudotimeKernel")
print("  - RealTimeKernel")
print("  - ConnectivityKernel")
print("  - GPCCA estimator")
```

---

## R Environment (Optional: Monocle3, Slingshot)

R tools are used in Week 2 for comparison with Python methods.

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

### Built-in Datasets (No Download Required)

```python
import scvelo as scv

# Pancreas — primary dataset for VelocityKernel labs
adata = scv.datasets.pancreas()

# Bone marrow — used for CytoTRACEKernel and PseudotimeKernel labs
adata = scv.datasets.bonemarrow()

# Dentate gyrus — alternative dataset
adata = scv.datasets.dentategyrus()
```

### CellRank Protocol Datasets (Figshare)

The CellRank Protocol provides curated datasets for all use cases:

```
https://doi.org/10.6084/m9.figshare.c.7752290.v1
```

These include preprocessed data for:
- CytoTRACEKernel (bone marrow hematopoiesis)
- PseudotimeKernel (bone marrow with DPT)
- RealTimeKernel (time-resolved data)
- VelocityKernel (pancreas, with metabolic labeling variant)

---

## Directory Structure

```
5.trajectory_analysis_course/
├── data/
│   ├── raw/
│   ├── velocity/
│   └── processed/
├── labs/
│   ├── lab03_diffusion_pt.ipynb
│   ├── lab06_velocity.ipynb
│   ├── lab07_trajectory_de.ipynb
│   ├── lab09_cellrank.ipynb          # VelocityKernel (Lab 8)
│   ├── lab09A_cellrank_cytotrace_kernel.ipynb  # CytoTRACEKernel (Lab 9A)
│   ├── lab09B_cellrank_pseudotime_kernel.ipynb # PseudotimeKernel (Lab 9B)
│   ├── lab09C_cellrank_realtime_kernel.ipynb   # RealTimeKernel (Lab 9C)
│   ├── lab09D_cellrank_kernel_combination.ipynb # Kernel Combination (Lab 9D)
│   └── lab10_comparison.ipynb
├── assignments/
├── quizzes/
├── results/
│   ├── pseudotime/
│   ├── velocity/
│   └── cellrank/
├── figures/
├── resources/
│   └── cellrank2_multiview_guide.md
└── cellrank_protocol/  ← cloned reference repo
```

---

## Computational Requirements

| Analysis | RAM | Time (5k cells) | GPU |
|----------|-----|-----------------|-----|
| Diffusion PT | 8 GB | 2 min | No |
| Monocle3 | 8 GB | 5 min | No |
| Slingshot | 4 GB | 2 min | No |
| scVelo stochastic | 8 GB | 5 min | No |
| scVelo dynamical | 16 GB | 30 min | No |
| CellRank VelocityKernel | 16 GB | 10 min | No |
| CellRank CytoTRACEKernel | 8 GB | 5 min | No |
| CellRank PseudotimeKernel | 8 GB | 5 min | No |
| CellRank RealTimeKernel | 16 GB | 15 min | No |
| CellRank kernel combination | 16 GB | 15 min | No |
| Full multiview pipeline | 16 GB | 60 min | No |

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `ModuleNotFoundError: cellrank` | `conda install -c conda-forge cellrank` |
| `ImportError: CytoTRACEKernel` | Update CellRank: `pip install --upgrade cellrank` (need v2.0+) |
| No velocity data | Need to process with velocyto/kb-python first |
| Monocle3 install fails | Check R version compatibility (≥4.1) |
| CellRank memory error | Use fewer PCs, fewer neighbors, or subsample |
| Velocity arrows chaotic | Try dynamical model; check spliced/unspliced quality |
| No clear trajectory | Data may not have a continuous process |
| GPCCA fails | Try fewer `n_components` in `compute_schur()` |
| Kernel combination error | Ensure all kernels use the same `adata` object |
| RealTimeKernel slow | Reduce `threshold` or subsample time points |

---

## Key References for Setup

- CellRank documentation: [cellrank.readthedocs.io](https://cellrank.readthedocs.io/)
- CellRank Protocol: [github.com/theislab/cellrank_protocol](https://github.com/theislab/cellrank_protocol)
- scVelo documentation: [scvelo.readthedocs.io](https://scvelo.readthedocs.io/)
- scanpy documentation: [scanpy.readthedocs.io](https://scanpy.readthedocs.io/)
