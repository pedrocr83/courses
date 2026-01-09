# Course: Clustering & Cell Type Annotation

## From Processed Counts to Biological Identity

---

## Course Overview

| | |
|---|---|
| **Duration** | 3 weeks (15-20 hours total) |
| **Format** | Self-paced with hands-on labs |
| **Level** | Beginner to Intermediate |
| **Prerequisites** | scRNA-seq processing, basic R/Python |

## Start Here (Do This First)

- **Step-by-step checklist**: `START_HERE.md`
- **Glossary**: `glossary.md`
- **Environment setup**: `setup.md`

### Learning Objectives

By the end of this course, you will be able to:
1. Apply dimensionality reduction techniques (PCA, UMAP, t-SNE)
2. Cluster cells using graph-based and other methods
3. Identify marker genes for each cluster
4. Annotate cell types using manual and automated approaches
5. Evaluate and refine annotations for downstream analysis

---

## Course Structure

### Week 1: Dimensionality Reduction

#### Module 1: Why Reduce Dimensions? (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 1.1 | The curse of dimensionality | Lecture |
| 1.2 | Goals: visualization vs analysis | Lecture |
| 1.3 | Overview of DR methods | Lecture |
| **Lab 1** | Explore high-dimensional gene expression data | Hands-on |

**Key Concepts:** Feature space, manifold hypothesis, information preservation

---

#### Module 2: Principal Component Analysis (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 2.1 | Linear algebra foundations | Lecture |
| 2.2 | PCA step-by-step | Lecture |
| 2.3 | Choosing number of PCs | Lecture |
| 2.4 | Highly variable genes selection | Lecture |
| **Lab 2** | Run PCA, interpret loadings, select PCs | Hands-on |

**Key Concepts:** Variance explained, elbow plot, HVGs, loadings

---

#### Module 3: Non-linear Methods (UMAP, t-SNE) (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 3.1 | t-SNE: perplexity and local structure | Lecture |
| 3.2 | UMAP: preserving global structure | Lecture |
| 3.3 | Comparing t-SNE vs UMAP | Lecture |
| 3.4 | Parameter tuning and interpretation | Demo |
| **Lab 3** | Generate and compare UMAP/t-SNE embeddings | Hands-on |

**Key Concepts:** Perplexity, n_neighbors, min_dist, reproducibility

---

### Week 2: Clustering

#### Module 4: Clustering Fundamentals (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 4.1 | What is a cluster? | Lecture |
| 4.2 | Distance metrics for gene expression | Lecture |
| 4.3 | Clustering approaches overview | Lecture |
| **Lab 4** | Explore cell-cell distances | Hands-on |

**Key Concepts:** Euclidean vs correlation distance, similarity graphs

---

#### Module 5: Graph-Based Clustering (2.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 5.1 | k-NN graphs and SNN graphs | Lecture |
| 5.2 | Louvain algorithm | Lecture |
| 5.3 | Leiden algorithm (improved) | Lecture |
| 5.4 | Resolution parameter and over-clustering | Lecture |
| **Lab 5A** | Build neighbor graphs and cluster | Hands-on |
| **Lab 5B** | Explore resolution effects | Hands-on |

**Key Concepts:** Community detection, resolution, modularity

---

#### Module 6: Evaluating Clusters (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 6.1 | Cluster quality metrics | Lecture |
| 6.2 | Silhouette scores and stability | Lecture |
| 6.3 | Biological validation | Discussion |
| **Lab 6** | Assess cluster quality and stability | Hands-on |

**Key Concepts:** Silhouette, ARI, marker separation

---

### Week 3: Cell Type Annotation

#### Module 7: Marker Gene Identification (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 7.1 | Differential expression for markers | Lecture |
| 7.2 | Statistical tests: Wilcoxon, t-test, logistic regression | Lecture |
| 7.3 | Effect size and specificity | Lecture |
| 7.4 | Visualizing markers: dotplots, violin, heatmaps | Demo |
| **Lab 7** | Find and visualize cluster markers | Hands-on |

**Key Concepts:** log2FC, pct.1/pct.2, specificity scores

---

#### Module 8: Manual Annotation (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 8.1 | Using known marker genes | Lecture |
| 8.2 | Literature and database resources | Lecture |
| 8.3 | Annotation workflow best practices | Demo |
| 8.4 | Handling ambiguous clusters | Discussion |
| **Lab 8** | Manually annotate PBMC clusters | Hands-on |

**Key Concepts:** Canonical markers, PanglaoDB, CellMarker, iterative refinement

---

#### Module 9: Automated Annotation (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 9.1 | Reference-based methods: SingleR, scmap | Lecture |
| 9.2 | Marker-based methods: CellTypist, scType | Lecture |
| 9.3 | Comparing and combining approaches | Lecture |
| **Lab 9** | Run automated annotation, compare to manual | Hands-on |

**Key Concepts:** Reference atlases, transfer learning, confidence scores

---

#### Module 10: Annotation Refinement & QC (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 10.1 | Sub-clustering for heterogeneity | Lecture |
| 10.2 | Merging over-split clusters | Lecture |
| 10.3 | Doublet contamination in clusters | Lecture |
| 10.4 | Documentation and reproducibility | Lecture |
| **Lab 10** | Refine annotations and create final labels | Hands-on |

**Key Concepts:** Hierarchical annotation, confidence thresholds, metadata

---

## Assessments

### Quizzes (after each module)
- 5-10 questions covering key concepts
- Immediate feedback

### Practical Assignments

| Assignment | Module | Description |
|------------|--------|-------------|
| **A1** | 2-3 | Dimensionality reduction comparison |
| **A2** | 5-6 | Clustering at multiple resolutions |
| **A3** | 7-9 | Complete cell type annotation |

### Final Project
**Annotate a novel scRNA-seq dataset end-to-end**
- Dimensionality reduction
- Clustering optimization
- Marker identification
- Manual + automated annotation
- Quality assessment report

---

## Software Setup

### Python (Scanpy)

```python
pip install scanpy anndata numpy pandas matplotlib seaborn
pip install leidenalg  # For Leiden clustering
pip install celltypist  # Automated annotation
```

### R (Seurat)

```r
install.packages("Seurat")
BiocManager::install(c("SingleR", "celldex", "scmap"))
install.packages(c("ggplot2", "pheatmap", "clustree"))
```

---

## Demo Datasets

| Dataset | Cells | Cell Types | Use |
|---------|-------|------------|-----|
| PBMC 3k | 2,700 | ~8 | Primary labs |
| PBMC 10k | 10,000 | ~12 | Assignment 3 |
| Tabula Muris | 100k+ | Many | Reference |

---

## Resources

### Marker Databases
| Resource | URL |
|----------|-----|
| PanglaoDB | https://panglaodb.se/ |
| CellMarker | http://bio-bigdata.hrbmu.edu.cn/CellMarker/ |
| Azimuth | https://azimuth.hubmapconsortium.org/ |

### Reference Atlases
- Human Cell Atlas
- Tabula Sapiens
- Mouse Cell Atlas

---

## Common Mistakes to Avoid

| Mistake | Consequence | Prevention |
|---------|-------------|------------|
| Using raw counts for DR | Dominated by library size | Normalize first |
| Too few HVGs | Miss rare cell types | Use 2000-5000 HVGs |
| Single resolution | Over/under-clustering | Test multiple resolutions |
| Ignoring doublet clusters | False cell types | Check doublet markers |
| Over-relying on UMAP | Misinterpret distances | Use for visualization only |
| No marker validation | Wrong annotations | Check known markers |

---

## What Bad Clustering / Bad Annotation Looks Like (and How to Catch It)

### Symptoms (Red Flags)

- **Embeddings**: UMAP changes drastically run-to-run; clusters appear/disappear with small parameter changes.
- **Clusters**: clusters separate by technical covariates (UMIs/mt%); clusters have no clear markers.
- **Markers**: “markers” are mostly mitochondrial/ribosomal; top markers are lowly expressed everywhere.
- **Annotation**: one cluster expresses canonical markers from multiple lineages; automated labels have low confidence.

### Common Failure Modes → What It Looks Like → Fix

- **Technical clustering (library size / mt%)**  
  - **Looks like**: clusters correlate with `total_counts`, `pct_counts_mt`; markers are QC genes.  
  - **Fix**: re-check normalization; regress technical covariates (carefully); remove low-quality cells; avoid using QC genes as HVGs.

- **Over-clustering (resolution too high)**  
  - **Looks like**: many tiny clusters; markers weak; “subtypes” not reproducible.  
  - **Fix**: lower resolution; require biological justification; use stability checks across resolutions.

- **Under-clustering (resolution too low)**  
  - **Looks like**: known cell types merged; mixed marker expression; broad clusters.  
  - **Fix**: increase resolution; subcluster specific populations; check marker gradients.

- **Doublet-driven clusters**  
  - **Looks like**: co-expression of lineage-incompatible markers (e.g., MS4A1 + NKG7).  
  - **Fix**: run doublet detection upstream; remove suspected doublets; re-cluster.

- **Annotation over-confidence**  
  - **Looks like**: assigning exact subtype names with weak evidence.  
  - **Fix**: label at the right granularity; record confidence; include marker evidence table.

### Minimum “Sanity Checks” Before You Proceed

- **Markers**: every cluster has at least a few plausible, specific markers.
- **Confounders**: clusters are not explained primarily by QC metrics.
- **Reproducibility**: clustering structure is stable across reasonable parameter ranges.
- **Documentation**: final label table includes evidence + confidence.

---

## Learning Path Progression

```
BEGINNER (Weeks 1-2)
│
├── Run PCA on scRNA-seq data
├── Generate UMAP embeddings
├── Cluster with Leiden
└── Find marker genes
        │
        ▼
INTERMEDIATE (Week 3)
│
├── Manual annotation with markers
├── Use automated tools (SingleR)
├── Assess annotation confidence
└── Sub-cluster heterogeneous populations
        │
        ▼
ADVANCED (Post-course)
│
├── Custom reference building
├── Cross-species annotation
├── Novel cell type discovery
└── Atlas-level integration
```

---

## Connection to Other Courses

| Prerequisite from | Why |
|-------------------|-----|
| scRNA-seq Processing | Need QC'd count matrix |

| Prerequisite for | Why |
|------------------|-----|
| Differential Expression | Compare annotated cell types |
| Cell-Cell Communication | CCC requires annotations |
| Trajectory Analysis | Annotated starting/end points |
| Data Integration | Need per-dataset annotations |

---

*This course is the critical bridge between data processing and biological interpretation.*

