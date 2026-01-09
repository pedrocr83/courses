# Assignment 4: Pseudobulk vs Cell-Level DE Comparison

**Modules:** 11-12 - Single-Cell DE Methods  
**Due:** End of Week 4  
**Points:** 150

---

## Overview

Compare pseudobulk and cell-level DE approaches on a single-cell dataset and evaluate the differences.

---

## Dataset

Use the PBMC dataset with simulated sample/donor information.

```python
import scanpy as sc
import numpy as np
import pandas as pd

# Load PBMC data
adata = sc.datasets.pbmc3k_processed()

# Simulate donors (biological replicates) for demonstration
np.random.seed(42)
n_cells = adata.n_obs

# Assign 6 donors (3 per condition for this exercise)
adata.obs['donor'] = np.random.choice(['D1', 'D2', 'D3', 'D4', 'D5', 'D6'], n_cells)
adata.obs['condition'] = adata.obs['donor'].map({
    'D1': 'control', 'D2': 'control', 'D3': 'control',
    'D4': 'treated', 'D5': 'treated', 'D6': 'treated'
})
```

*Note: This uses simulated conditions for educational purposes.*

---

## Part A: Cell-Level DE (30 points)

### Task 1: Wilcoxon Test

Using Scanpy, run cell-level Wilcoxon tests comparing conditions.

```python
sc.tl.rank_genes_groups(adata, groupby='condition', method='wilcoxon')
```

**Report:**
1. Number of DEGs (FDR < 0.05)
2. Number of DEGs (FDR < 0.05, |log2FC| > 0.5)
3. Top 20 DEGs

### Task 2: Critique

Discuss (3-4 sentences):
- What statistical assumption is violated?
- Why might results be unreliable?

---

## Part B: Pseudobulk Aggregation (40 points)

### Task 3: Aggregate Counts

Aggregate counts by donor and cell type (choose one cell type for simplicity).

```python
# Filter to one cell type (e.g., CD4 T cells)
adata_subset = adata[adata.obs['louvain'] == '0'].copy()

# Aggregate by donor
pseudobulk = adata_subset.to_df().groupby(adata_subset.obs['donor']).sum()
```

### Task 4: Create Metadata

Build sample metadata for the pseudobulk matrix.

```python
sample_meta = pd.DataFrame({
    'donor': pseudobulk.index,
    'condition': ['control', 'control', 'control', 'treated', 'treated', 'treated']
})
```

### Task 5: Run DESeq2-style Analysis

Option A (Python with pyDESeq2):
```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
```

Option B (Export to R):
```python
pseudobulk.T.to_csv('pseudobulk_counts.csv')
sample_meta.to_csv('sample_metadata.csv')
```

Then in R:
```r
library(DESeq2)
cts <- read.csv('pseudobulk_counts.csv', row.names=1)
coldata <- read.csv('sample_metadata.csv', row.names=1)
dds <- DESeqDataSetFromMatrix(cts, coldata, ~ condition)
dds <- DESeq(dds)
res <- results(dds)
```

---

## Part C: Comparison (40 points)

### Task 6: Compare Results

Create a comparison of the two approaches:

1. **Venn diagram** of DEGs from each method
2. **Scatter plot** of log2FC: cell-level vs pseudobulk
3. **Scatter plot** of -log10(p-value): cell-level vs pseudobulk

### Task 7: Quantify Differences

| Metric | Cell-Level | Pseudobulk |
|--------|------------|------------|
| Total DEGs (FDR < 0.05) | | |
| DEGs unique to method | | |
| Overlap | | |
| Correlation of log2FC | | |

---

## Part D: Discussion (40 points)

Write 400-600 words addressing:

1. **Which method identified more DEGs?** Why?

2. **Which results are more trustworthy?** Why?

3. **When would cell-level tests be acceptable?**
   - Exploratory analysis?
   - Marker gene finding?
   - Formal hypothesis testing?

4. **What are the limitations of pseudobulk?**
   - Cell type heterogeneity
   - Minimum cells per sample
   - Loss of cell-level resolution

5. **Recommendations:**
   - What would you recommend for a collaborator's scRNA-seq DE analysis?

---

## Deliverables

1. **Jupyter notebook or R Markdown** with all code
2. **Figures:**
   - Venn diagram
   - log2FC comparison scatter
   - p-value comparison scatter
3. **Written discussion** (separate section or document)
4. **CSV files:**
   - Cell-level DE results
   - Pseudobulk DE results

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Cell-level analysis correct | 30 |
| Pseudobulk aggregation correct | 30 |
| DESeq2/pyDESeq2 run correctly | 20 |
| Comparison visualizations | 30 |
| Discussion depth | 30 |
| Code quality | 10 |

---

## Key Takeaways

By completing this assignment, you should understand:
- Cell-level tests inflate significance due to treating cells as independent
- Pseudobulk respects biological replication
- The number of DEGs is not a measure of method quality
- Appropriate method depends on experimental design and goals

