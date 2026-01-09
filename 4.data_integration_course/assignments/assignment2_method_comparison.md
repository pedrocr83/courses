# Assignment 2: Harmony vs Seurat Integration

**Modules:** 5-6 - Harmony & Seurat  
**Due:** End of Week 2  
**Points:** 150

---

## Overview

Compare two popular integration methods on the same dataset.

---

## Tasks

### Part A: Run Harmony (40 points)

1. Run Harmony integration
2. Document parameters used
3. Generate UMAP from corrected embedding
4. Cluster the integrated data
5. Visualize by batch and cell type

---

### Part B: Run Seurat Integration (40 points)

1. Run Seurat CCA or RPCA integration
2. Document parameters and anchor selection
3. Generate UMAP
4. Cluster the integrated data
5. Visualize by batch and cell type

---

### Part C: Quantitative Comparison (40 points)

For both methods, calculate:
- iLISI (batch mixing)
- cLISI (cell type conservation)
- ARI (cluster vs known labels)

Create comparison table:

| Metric | Pre-integration | Harmony | Seurat |
|--------|-----------------|---------|--------|
| iLISI | | | |
| cLISI | | | |
| ARI | | | |

---

### Part D: Discussion (30 points)

Write 300 words:

1. Which method achieved better batch mixing?
2. Which preserved biology better?
3. Which would you recommend for this dataset and why?
4. What are the trade-offs between methods?

---

## Deliverables

1. Notebooks for both methods
2. Figures: UMAPs before/after for each method
3. Metrics comparison table
4. Written discussion

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Harmony analysis | 40 |
| Seurat analysis | 40 |
| Quantitative comparison | 40 |
| Discussion | 30 |

