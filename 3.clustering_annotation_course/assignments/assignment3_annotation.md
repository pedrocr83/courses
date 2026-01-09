# Assignment 3: Complete Cell Type Annotation

**Modules:** 7-9 - Marker Genes & Annotation  
**Due:** End of Week 3  
**Points:** 150

---

## Overview

Perform complete cell type annotation using both manual and automated approaches.

---

## Tasks

### Part A: Find Marker Genes (30 points)

1. Run differential expression to find markers for each cluster
2. For each cluster, identify top 5 marker genes
3. Create:
   - Dot plot of top 3 markers per cluster
   - Heatmap of top markers

---

### Part B: Manual Annotation (40 points)

Using known PBMC markers:

| Cell Type | Key Markers |
|-----------|-------------|
| CD4 T cells | CD3D, CD4, IL7R |
| CD8 T cells | CD3D, CD8A, CD8B |
| B cells | CD19, MS4A1, CD79A |
| NK cells | NKG7, GNLY, NCAM1 |
| CD14 Monocytes | CD14, LYZ, S100A8 |
| CD16 Monocytes | FCGR3A, MS4A7 |
| Dendritic cells | FCER1A, CST3 |
| Platelets | PPBP, PF4 |

1. Check expression of these markers in your clusters
2. Assign cell type labels to each cluster
3. Document any ambiguous clusters and your reasoning

---

### Part C: Automated Annotation (30 points)

1. Run SingleR or CellTypist on your data
2. Compare automated labels to your manual annotations
3. Calculate agreement percentage

---

### Part D: Final Annotation (30 points)

1. Create final cell type assignments
2. Where manual and automated disagree, investigate and decide
3. Create publication-quality UMAP with cell type labels
4. Create summary table:

| Cluster | Manual | Automated | Final | Cells |
|---------|--------|-----------|-------|-------|
| 0 | | | | |

---

### Part E: Reflection (20 points)

Write 200-300 words:

1. Where did manual and automated annotation agree/disagree?
2. What were the most challenging clusters to annotate?
3. How confident are you in each annotation?

---

## Deliverables

1. Notebook with full analysis
2. Figures: Marker dot plot, heatmap, final UMAP
3. Annotation summary table
4. Written reflection

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Marker identification | 30 |
| Manual annotation | 40 |
| Automated annotation | 30 |
| Final annotation | 30 |
| Reflection quality | 20 |

