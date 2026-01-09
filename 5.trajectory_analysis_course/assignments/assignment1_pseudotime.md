# Assignment 1: Pseudotime Inference Comparison

**Modules:** 3-4 - Diffusion PT & Monocle  
**Due:** End of Week 1  
**Points:** 100

---

## Overview

Compare pseudotime inference methods on developmental data.

---

## Dataset

Use the pancreas endocrine differentiation dataset.

```python
import scvelo as scv
adata = scv.datasets.pancreas()
```

---

## Tasks

### Part A: Data Exploration (20 points)

1. Explore the dataset structure
2. Identify annotated cell types
3. Visualize on UMAP
4. Identify expected developmental progression

---

### Part B: Diffusion Pseudotime (30 points)

1. Compute diffusion map
2. Identify root cell (use biological knowledge)
3. Compute diffusion pseudotime
4. Visualize pseudotime on UMAP
5. Check if ordering matches expected biology

---

### Part C: Monocle 3 (30 points)

1. Create Monocle3 object
2. Learn principal graph
3. Order cells (choose root)
4. Visualize trajectory
5. Compare to diffusion PT

---

### Part D: Comparison (20 points)

1. Correlate pseudotime values between methods
2. Identify agreements and disagreements
3. Write 150 words on which method better captures the biology

---

## Deliverables

1. Notebook with both methods
2. Figures: UMAP colored by pseudotime (both methods)
3. Correlation plot
4. Written comparison

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Data exploration | 20 |
| Diffusion PT | 30 |
| Monocle 3 | 30 |
| Comparison | 20 |

