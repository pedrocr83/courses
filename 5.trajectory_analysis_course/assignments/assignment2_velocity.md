# Assignment 2: RNA Velocity Analysis

**Module:** 6 - RNA Velocity  
**Due:** End of Week 2  
**Points:** 100

---

## Overview

Compute and interpret RNA velocity to understand differentiation direction.

---

## Dataset

Use a dataset with velocity information (spliced/unspliced counts).

```python
import scvelo as scv
adata = scv.datasets.pancreas()
# or
adata = scv.datasets.dentategyrus()
```

---

## Tasks

### Part A: Velocity Computation (30 points)

1. Check spliced/unspliced layers exist
2. Filter and normalize for velocity
3. Compute moments
4. Run velocity (stochastic model)

```python
scv.pp.filter_and_normalize(adata, min_shared_counts=20)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
```

---

### Part B: Visualization (30 points)

1. Create velocity stream plot on UMAP
2. Create velocity arrows on UMAP
3. Plot velocity for specific genes
4. Identify regions of high vs low velocity

---

### Part C: Interpretation (20 points)

1. Which cell types show highest velocity?
2. What directions do arrows point?
3. Does this match expected differentiation?
4. Identify any unexpected patterns

---

### Part D: Advanced (20 points)

1. Run dynamical model and compare
2. Identify velocity genes (high vs low R²)
3. Discuss limitations of velocity in this dataset

---

## Deliverables

1. Notebook with velocity analysis
2. Figures: Stream plot, arrow plot, gene plots
3. Written interpretation (200 words)

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Velocity computation | 30 |
| Visualization | 30 |
| Interpretation | 20 |
| Advanced analysis | 20 |

