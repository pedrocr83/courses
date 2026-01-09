# Assignment 1: Prepare Data for CCC Analysis

**Modules:** 2-3 - Data Preparation  
**Due:** End of Week 1  
**Points:** 100

---

## Overview

Prepare a scRNA-seq dataset for cell-cell communication analysis, ensuring proper annotation and formatting.

---

## Dataset

Use the PBMC 3k dataset (or another approved dataset).

```python
import scanpy as sc
adata = sc.datasets.pbmc3k_processed()
```

---

## Tasks

### Part A: Data Inspection (25 points)

1. Load the dataset and describe its properties:
   - Number of cells
   - Number of genes
   - Cell type annotations present
   - Conditions/batches if any

2. Visualize cell type composition:
   - Bar plot of cell type frequencies
   - UMAP colored by cell type

3. Check annotation quality:
   - Are cell types well-separated on UMAP?
   - Are marker genes expressed appropriately?

---

### Part B: Expression Thresholds (25 points)

1. Calculate what percentage of genes are expressed (count > 0) in each cell type
2. Plot distribution of mean expression per cell type
3. Recommend a minimum expression threshold for L-R analysis
4. Justify your threshold choice

---

### Part C: Format for CCC Tools (25 points)

1. **For CellPhoneDB (Python):**
   - Export counts matrix (genes × cells)
   - Export metadata with cell type annotations
   - Verify format requirements

2. **For CellChat (R):**
   - Convert to Seurat object (if needed)
   - Ensure proper assay structure
   - Set default identity to cell types

```python
# CellPhoneDB format
adata.obs[['cell_type']].to_csv('meta.txt', sep='\t')
adata.to_df().T.to_csv('counts.txt', sep='\t')
```

```r
# CellChat format
library(CellChat)
cellchat <- createCellChat(object = seurat_obj, group.by = "cell_type")
```

---

### Part D: Ligand-Receptor Check (25 points)

1. Pick 5 known L-R pairs relevant to immune cells:
   - e.g., CD40LG-CD40, CCL5-CCR5, IL6-IL6R

2. Check their expression in your data:
   - Which cell types express the ligand?
   - Which cell types express the receptor?
   - Is the expression above your threshold?

3. Create a simple table of expected interactions based on expression

| Ligand | Ligand+ Cell Types | Receptor | Receptor+ Cell Types |
|--------|-------------------|----------|---------------------|
| | | | |

---

## Deliverables

1. **Jupyter notebook or R script** with all code
2. **Figures:**
   - UMAP with annotations
   - Cell type bar plot
   - Expression distribution
3. **Exported files:**
   - CellPhoneDB-ready files
   - CellChat-ready object (or code to create it)
4. **L-R expression table**

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Data inspection complete | 25 |
| Threshold analysis thoughtful | 25 |
| Correct file formats | 25 |
| L-R expression check | 25 |

---

## Tips

- CCC results are only as good as your annotations
- Rare cell types may not have enough power
- Expression threshold affects sensitivity vs specificity

