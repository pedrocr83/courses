# Assignment 1: Dimensionality Reduction Comparison

**Modules:** 2-3 - PCA, UMAP, t-SNE  
**Due:** End of Week 1  
**Points:** 100

---

## Overview

Compare dimensionality reduction methods on scRNA-seq data and understand their strengths and limitations.

---

## Dataset

Use the PBMC 3k dataset.

```python
import scanpy as sc
adata = sc.datasets.pbmc3k()
```

---

## Tasks

### Part A: Data Preparation (20 points)

1. Load and perform basic QC filtering
2. Normalize (log1p)
3. Select highly variable genes (2000 HVGs)
4. Document number of cells and genes after filtering

---

### Part B: PCA Analysis (30 points)

1. Run PCA
2. Create elbow plot (variance explained)
3. Determine optimal number of PCs
4. Examine PC1 and PC2 loadings - which genes contribute most?
5. Check if PC1 correlates with total counts (technical variation)

---

### Part C: UMAP vs t-SNE (30 points)

1. Run UMAP with default parameters
2. Run t-SNE with default parameters
3. Run UMAP with different `n_neighbors` (5, 15, 50)
4. Run t-SNE with different `perplexity` (5, 30, 100)

Create a figure grid showing all embeddings.

---

### Part D: Critical Comparison (20 points)

Write 200-300 words addressing:

1. How do UMAP and t-SNE differ visually?
2. How do parameters affect cluster appearance?
3. What can and cannot be inferred from these plots?
4. Which would you use for publication and why?

---

## Deliverables

1. Jupyter notebook with all code
2. Figure: 6-panel comparison of embeddings
3. Written comparison

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Data preparation correct | 20 |
| PCA analysis complete | 30 |
| Multiple embeddings generated | 30 |
| Thoughtful comparison | 20 |

