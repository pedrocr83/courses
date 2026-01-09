# Assignment 3: NicheNet Ligand Activity + LIANA Consensus

**Modules:** 6-7 - NicheNet & LIANA  
**Due:** End of Week 3  
**Points:** 150

---

## Overview

Use NicheNet to predict ligand activities and LIANA to generate consensus rankings across multiple methods.

---

## Dataset

Use PBMC data with defined sender and receiver cell types.

---

## Part A: NicheNet Analysis (60 points)

### Setup

```r
library(nichenetr)
library(Seurat)

# Load prior networks (download once)
ligand_target_matrix <- readRDS("ligand_target_matrix.rds")
lr_network <- readRDS("lr_network.rds")
```

### Task 1: Define the Biological Question (15 points)

1. Choose a sender cell type (e.g., T cells)
2. Choose a receiver cell type (e.g., Monocytes)
3. Define genes of interest in receiver (e.g., DEGs between conditions, or marker genes)

Document your choices and biological rationale.

### Task 2: Run NicheNet (25 points)

```r
# Define expressed genes
expressed_genes_sender <- get_expressed_genes(sender_celltype, seurat_obj)
expressed_genes_receiver <- get_expressed_genes(receiver_celltype, seurat_obj)

# Define genes of interest (target genes in receiver)
geneset_oi <- DEGs_receiver  # or marker genes

# Define background
background_genes <- expressed_genes_receiver

# Predict ligand activities
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)
```

### Task 3: Interpret Results (20 points)

1. Identify top 10 ligands by activity score (Pearson correlation)
2. For top 3 ligands, identify:
   - Target genes in receiver
   - Receptors in receiver
3. Visualize ligand-target connections

---

## Part B: LIANA Consensus (50 points)

### Task 4: Run LIANA with Multiple Methods (25 points)

```r
library(liana)

# Run LIANA with all methods
liana_results <- liana_wrap(seurat_obj)

# View methods used
liana_results
```

Or in Python:

```python
import liana as li

# Run LIANA
li.mt.cellphonedb(adata)
li.mt.rank_aggregate(adata)
```

### Task 5: Analyze Consensus (25 points)

1. Extract aggregate rankings
2. Identify top 20 consensus interactions
3. Compare to individual method rankings
4. Create heatmap of method agreement

---

## Part C: Integration (40 points)

### Task 6: Compare NicheNet and LIANA Results (20 points)

1. Do LIANA's top interactions involve NicheNet's top ligands?
2. Create a table showing:
   - NicheNet ligand rank
   - LIANA interaction rank (if that ligand is involved)
   - Agreement/disagreement

### Task 7: Biological Summary (20 points)

Write 200-300 words:
1. Which ligands are consistently identified?
2. What biological processes do they regulate?
3. How do NicheNet's target predictions add to L-R inference?

---

## Deliverables

1. **R notebook** with NicheNet and LIANA code
2. **Figures:**
   - NicheNet ligand activity barplot
   - NicheNet ligand-target heatmap
   - LIANA consensus heatmap
   - Method agreement visualization
3. **Tables:**
   - Top ligands (NicheNet)
   - Top interactions (LIANA)
   - Comparison table
4. **Written interpretation**

---

## Rubric

| Criterion | Points |
|-----------|--------|
| NicheNet setup and rationale | 15 |
| NicheNet analysis complete | 25 |
| NicheNet interpretation | 20 |
| LIANA multi-method run | 25 |
| LIANA consensus analysis | 25 |
| Integration of both tools | 20 |
| Biological interpretation | 20 |

---

## Key Takeaways

- NicheNet adds ligand → target causality (not just co-expression)
- LIANA reduces method-specific bias
- Combining approaches increases confidence

