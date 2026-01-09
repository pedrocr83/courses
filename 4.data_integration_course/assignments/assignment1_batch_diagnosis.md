# Assignment 1: Diagnose and Visualize Batch Effects

**Modules:** 2-3 - Diagnosing Batch Effects  
**Due:** End of Week 1  
**Points:** 100

---

## Overview

Learn to identify and quantify batch effects in multi-sample scRNA-seq data.

---

## Dataset

Use the pancreas integration benchmark dataset (4 technologies).

---

## Tasks

### Part A: Load and Explore (20 points)

1. Load the multi-batch dataset
2. Document:
   - Number of cells per batch
   - Number of batches
   - Cell type annotations (if available)
3. Process each batch identically (normalize, HVG, PCA)

---

### Part B: Visual Diagnosis (30 points)

1. Create combined UMAP colored by:
   - Batch/technology
   - Cell type (if available)
2. Describe what you observe:
   - Do batches cluster separately?
   - Are cell types separated across batches?

---

### Part C: Quantitative Metrics (30 points)

1. Calculate LISI scores (iLISI and cLISI)
2. Calculate kBET acceptance rate (if available)
3. Interpret the metrics:
   - What do the values tell you about mixing?
   - What do they tell you about biology preservation?

---

### Part D: Assessment (20 points)

Write 200 words addressing:

1. Does this data need integration?
2. What batch effects are present?
3. Is there any confounding between batch and biology?

---

## Deliverables

1. Notebook with analysis
2. Figure: UMAP by batch and cell type
3. Metrics table
4. Written assessment

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Data exploration | 20 |
| Visual diagnosis | 30 |
| Quantitative metrics | 30 |
| Assessment quality | 20 |

