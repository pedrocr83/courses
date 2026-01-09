# Assignment 2: Clustering at Multiple Resolutions

**Modules:** 5-6 - Graph-Based Clustering  
**Due:** End of Week 2  
**Points:** 100

---

## Overview

Explore how resolution affects clustering and learn to choose appropriate granularity.

---

## Tasks

### Part A: Build Neighbor Graph (20 points)

1. Load preprocessed PBMC data (from Assignment 1)
2. Compute k-NN graph with k=15
3. Verify graph was constructed correctly

---

### Part B: Multi-Resolution Clustering (40 points)

Run Leiden clustering at resolutions: 0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0

For each resolution, record:
- Number of clusters
- Sizes of clusters

Create:
1. UMAP colored by clusters for each resolution
2. Summary table of cluster counts

---

### Part C: Cluster Stability (20 points)

1. Use clustree (R) or manual tracking to visualize how clusters split/merge
2. Identify which clusters are stable across resolutions
3. Identify which clusters only appear at high resolution

---

### Part D: Choose Optimal Resolution (20 points)

Based on your analysis:

1. Recommend an optimal resolution for this dataset
2. Justify your choice with:
   - Number of expected cell types in PBMCs
   - Cluster stability
   - Biological interpretability (check a few markers)

---

## Deliverables

1. Notebook with clustering analysis
2. Figure: Multi-panel UMAP by resolution
3. Clustree visualization or equivalent
4. Written justification (150-200 words)

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Graph construction | 20 |
| Multi-resolution clustering | 40 |
| Stability analysis | 20 |
| Resolution justification | 20 |

