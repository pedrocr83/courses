# Course: Data Integration & Batch Correction

## Combining Multiple Datasets into Unified Atlases

---

## Course Overview

| | |
|---|---|
| **Duration** | 3 weeks (15-20 hours total) |
| **Format** | Self-paced with hands-on labs |
| **Level** | Intermediate |
| **Prerequisites** | scRNA-seq processing, clustering & annotation |

## Start Here (Do This First)

- **Step-by-step checklist**: `START_HERE.md`
- **Glossary**: `glossary.md`
- **Environment setup**: `setup.md`

### Learning Objectives

By the end of this course, you will be able to:
1. Identify and diagnose batch effects in scRNA-seq data
2. Apply linear and non-linear batch correction methods
3. Integrate datasets across conditions, technologies, and species
4. Evaluate integration quality with quantitative metrics
5. Build multi-sample atlases for downstream analysis

---

## Course Structure

### Week 1: Understanding Batch Effects

#### Module 1: What Are Batch Effects? (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 1.1 | Technical vs biological variation | Lecture |
| 1.2 | Sources: sequencing, processing, reagents | Lecture |
| 1.3 | Impact on downstream analysis | Case studies |
| **Lab 1** | Visualize batch effects in multi-sample data | Hands-on |

**Key Concepts:** Confounding, technical artifacts, experimental design

---

#### Module 2: Diagnosing Batch Effects (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 2.1 | Visual inspection: UMAP by batch | Lecture |
| 2.2 | Quantitative metrics: LISI, kBET | Lecture |
| 2.3 | When is correction needed? | Discussion |
| **Lab 2** | Compute and interpret batch metrics | Hands-on |

**Key Concepts:** Local Inverse Simpson Index, batch mixing, over-correction

---

#### Module 3: Simple Correction Methods (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 3.1 | Regressing out batch | Lecture |
| 3.2 | ComBat for bulk/pseudo-bulk | Lecture |
| 3.3 | Limitations of linear methods | Lecture |
| **Lab 3** | Apply regression-based correction | Hands-on |

**Key Concepts:** Linear regression, covariates, residuals

---

### Week 2: Integration Methods

#### Module 4: Mutual Nearest Neighbors (MNN) (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 4.1 | MNN concept and algorithm | Lecture |
| 4.2 | fastMNN in R | Demo |
| 4.3 | When MNN works well | Discussion |
| **Lab 4** | Run MNN integration | Hands-on |

**Key Concepts:** Anchor cells, batch vectors, local correction

---

#### Module 5: Harmony (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 5.1 | Iterative clustering approach | Lecture |
| 5.2 | PCA embedding correction | Lecture |
| 5.3 | Speed and scalability | Lecture |
| **Lab 5** | Integrate with Harmony | Hands-on |

**Key Concepts:** Soft clustering, iterative adjustment, embedding space

---

#### Module 6: Seurat Integration (CCA/RPCA) (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 6.1 | Canonical Correlation Analysis | Lecture |
| 6.2 | Anchor identification | Lecture |
| 6.3 | CCA vs RPCA: when to use each | Lecture |
| **Lab 6** | Seurat v5 integration workflow | Hands-on |

**Key Concepts:** Anchors, integration features, reference mapping

---

#### Module 7: Deep Learning Methods (scVI, scANVI) (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 7.1 | Variational autoencoders for scRNA-seq | Lecture |
| 7.2 | scVI: unsupervised integration | Lecture |
| 7.3 | scANVI: semi-supervised with labels | Lecture |
| **Lab 7** | Run scVI integration | Hands-on |

**Key Concepts:** Latent space, VAE, batch as covariate

---

### Week 3: Evaluation & Applications

#### Module 8: Evaluating Integration Quality (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 8.1 | Batch mixing metrics (iLISI, kBET) | Lecture |
| 8.2 | Biological conservation (cLISI, ARI) | Lecture |
| 8.3 | scIB benchmarking framework | Lecture |
| **Lab 8** | Comprehensive integration evaluation | Hands-on |

**Key Concepts:** Mixing vs conservation trade-off, benchmark scores

---

#### Module 9: Choosing the Right Method (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 9.1 | Method comparison studies | Reading |
| 9.2 | Decision tree for method selection | Lecture |
| 9.3 | When integration fails | Discussion |
| **Lab 9** | Compare 3 methods on same dataset | Hands-on |

**Key Concepts:** No free lunch, dataset-specific performance

---

#### Module 10: Reference Mapping & Label Transfer (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 10.1 | Mapping queries to references | Lecture |
| 10.2 | Label transfer strategies | Lecture |
| 10.3 | Azimuth and automated mapping | Demo |
| **Lab 10** | Map new data to PBMC reference | Hands-on |

**Key Concepts:** Reference atlases, projection, transfer learning

---

## Assessments

### Quizzes (after each module)
- 5-10 questions covering key concepts
- Immediate feedback

### Practical Assignments

| Assignment | Module | Description |
|------------|--------|-------------|
| **A1** | 2-3 | Diagnose and visualize batch effects |
| **A2** | 5-6 | Compare Harmony vs Seurat integration |
| **A3** | 8 | Full integration evaluation report |

### Final Project
**Integrate a multi-batch dataset and evaluate**
- Multiple samples/conditions
- Apply 2+ integration methods
- Quantitative comparison
- Downstream analysis on integrated data

---

## Software Setup

### Python (scvi-tools, Scanpy)

```python
pip install scanpy anndata scvi-tools harmonypy
pip install scib  # For benchmarking
```

### R (Seurat, Harmony)

```r
install.packages("Seurat")
install.packages("harmony")
BiocManager::install("batchelor")  # For MNN
```

---

## Demo Datasets

| Dataset | Batches | Description | Use |
|---------|---------|-------------|-----|
| Pancreas (4 datasets) | 4 technologies | Classic benchmark | Labs 4-9 |
| PBMC multi-donor | 8 donors | Biological replicates | Assignment |
| Cross-species | Human + Mouse | Advanced | Final project option |

---

## Integration Methods Summary

| Method | Approach | Speed | Best For |
|--------|----------|-------|----------|
| Harmony | Iterative soft clustering | Fast | Large datasets |
| Seurat CCA | Anchor-based | Medium | Similar datasets |
| Seurat RPCA | Reference projection | Medium | Query → Reference |
| MNN/fastMNN | Mutual nearest neighbors | Medium | Few shared cell types |
| scVI | Deep learning VAE | Slow (GPU) | Complex batch structure |
| scANVI | Semi-supervised VAE | Slow (GPU) | When labels available |

---

## Common Mistakes to Avoid

| Mistake | Consequence | Prevention |
|---------|-------------|------------|
| Integrating before QC | Propagate bad cells | QC each batch first |
| Over-correcting | Lose biological signal | Check cell type separation |
| Ignoring confounding | False biological conclusions | Document batch-condition overlap |
| Single method | Miss best approach | Compare multiple methods |
| No evaluation | Unknown quality | Always compute metrics |
| Integrating incompatible data | Forced false similarities | Check overlap in cell types |

---

## What Bad Integration Looks Like (and How to Catch It)

### Symptoms (Red Flags)

- **Before integration**: UMAPs separate by batch even within the same cell type.
- **After integration**: cell types collapse into each other (“soup”); known markers lose specificity.
- **Confounding**: “perfect” mixing but condition differences vanish; batch and condition are indistinguishable.
- **Metrics**: improved batch mixing but destroyed biological conservation (trade-off ignored).

### Common Failure Modes → What It Looks Like → Fix

- **Over-correction**  
  - **Looks like**: distinct cell types blend; markers no longer separate populations.  
  - **Fix**: reduce correction strength; integrate within cell types; consider reference mapping; evaluate conservation metrics.

- **Under-correction**  
  - **Looks like**: batches remain separated within the same cell type.  
  - **Fix**: try different method; tune parameters; ensure embeddings built from comparable preprocessing.

- **Integrating incompatible datasets**  
  - **Looks like**: forced alignment of cell types that don’t exist in both datasets; spurious “matches”.  
  - **Fix**: integrate only shared compartments; remove unmatched populations or treat as query-only; use mapping not full integration.

- **Confounded design (batch=condition)**  
  - **Looks like**: cannot tell if differences are technical or biological; correction may erase real signal.  
  - **Fix**: redesign; if impossible, avoid strong claims; show sensitivity analyses.

### Minimum “Sanity Checks” Before Downstream Analysis

- **Within each cell type**: batches should mix; across cell types, biology should remain separated.
- **Markers**: canonical markers still define populations post-integration.
- **Metrics**: report at least one mixing metric + one conservation metric.
- **Documentation**: record method + parameters + preprocessing parity across batches.

---

## Learning Path Progression

```
BEGINNER (Week 1)
│
├── Visualize batch effects
├── Compute LISI metrics
└── Simple regression correction
        │
        ▼
INTERMEDIATE (Week 2)
│
├── Harmony integration
├── Seurat CCA/RPCA
├── Evaluate mixing vs conservation
└── Choose method for dataset
        │
        ▼
ADVANCED (Week 3+)
│
├── scVI/scANVI deep learning
├── Reference atlas mapping
├── Cross-species integration
└── Building new references
```

---

## Connection to Other Courses

| Prerequisite from | Why |
|-------------------|-----|
| scRNA-seq Processing | Need QC'd data |
| Clustering & Annotation | Need per-batch annotations |

| Prerequisite for | Why |
|------------------|-----|
| Differential Expression | Compare across integrated conditions |
| Cell-Cell Communication | Integrated cell types |
| Trajectory Analysis | Unified developmental space |

---

*This course enables analysis of real-world multi-sample experiments.*

