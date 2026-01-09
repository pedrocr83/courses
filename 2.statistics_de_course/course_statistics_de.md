# Course: Statistics & Differential Expression for Transcriptomics

## From Count Data to Biological Insight

---

## Course Overview

| | |
|---|---|
| **Duration** | 4 weeks (20-25 hours total) |
| **Format** | Self-paced with hands-on labs |
| **Level** | Beginner to Intermediate |
| **Prerequisites** | Basic R or Python, familiarity with gene expression concepts |

## Start Here (Do This First)

- **Step-by-step checklist**: `START_HERE.md`
- **Glossary**: `glossary.md`
- **Environment setup**: `setup.md`

### Learning Objectives

By the end of this course, you will be able to:
1. Understand probability distributions underlying RNA-seq data
2. Perform differential expression analysis on bulk RNA-seq
3. Apply pseudobulk DE methods for single-cell data
4. Correctly interpret DE results and avoid common pitfalls
5. Visualize and communicate findings effectively

---

## Course Structure

### Week 1: Foundations of Transcriptomics Statistics

#### Module 1: Introduction to Transcriptomics Data (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 1.1 | What is transcriptomics? | Lecture |
| 1.2 | Bulk vs single-cell vs spatial | Lecture |
| 1.3 | The count matrix: structure and properties | Lecture |
| **Lab 1** | Explore a count matrix in R/Python | Hands-on |

**Key Concepts:** RNA abundance, count data, mean-variance relationship

---

#### Module 2: Probability Distributions for Count Data (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 2.1 | Why we need probability models | Lecture |
| 2.2 | Poisson distribution: simple but limited | Lecture |
| 2.3 | Negative binomial: handling overdispersion | Lecture |
| 2.4 | Zero-inflated models for scRNA-seq | Reading |
| **Lab 2** | Visualize and fit distributions to RNA-seq data | Hands-on |

**Key Concepts:** Overdispersion, variance modeling, distribution fitting

---

#### Module 3: Experimental Design (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 3.1 | Biological vs technical replicates | Lecture |
| 3.2 | Confounders: batch, donor, protocol | Lecture |
| 3.3 | Design matrices explained | Lecture + Demo |
| 3.4 | Power analysis basics | Reading |
| **Lab 3** | Build design matrices for different experiments | Hands-on |

**Key Concepts:** Replication, confounding, balanced designs

---

### Week 2: Normalization & Hypothesis Testing

#### Module 4: Normalization Methods (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 4.1 | Why raw counts are misleading | Lecture |
| 4.2 | Library size normalization | Lecture |
| 4.3 | TMM and median-of-ratios | Lecture |
| 4.4 | scRNA-seq normalization: CP10K, SCTransform | Lecture |
| **Lab 4** | Compare normalization methods on real data | Hands-on |

**Key Concepts:** Sequencing depth, size factors, variance stabilization

---

#### Module 5: Hypothesis Testing Fundamentals (2.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 5.1 | Null hypothesis and test statistics | Lecture |
| 5.2 | p-values: what they mean and don't mean | Lecture |
| 5.3 | Type I and Type II errors | Lecture |
| 5.4 | Common statistical tests for expression | Lecture |
| **Lab 5** | Run t-tests and Wilcoxon tests on gene expression | Hands-on |

**Key Concepts:** H₀, p-value interpretation, test selection

---

#### Module 6: Multiple Testing Correction (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 6.1 | The multiple testing problem | Lecture |
| 6.2 | Bonferroni correction: too strict? | Lecture |
| 6.3 | False Discovery Rate (FDR) and Benjamini-Hochberg | Lecture |
| **Lab 6** | Apply FDR correction and compare results | Hands-on |

**Key Concepts:** FWER, FDR, adjusted p-values

---

### Week 3: Differential Expression Methods

#### Module 7: Effect Size & Biological Meaning (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 7.1 | Statistical vs biological significance | Lecture |
| 7.2 | Log2 fold change explained | Lecture |
| 7.3 | Shrinkage estimators | Lecture |
| **Lab 7** | Calculate and interpret effect sizes | Hands-on |

**Key Concepts:** log2FC, expression thresholds, shrinkage

---

#### Module 8: Bulk RNA-seq DE with DESeq2 (3 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 8.1 | DESeq2 workflow overview | Lecture |
| 8.2 | Negative binomial GLM | Lecture |
| 8.3 | Running DESeq2 step-by-step | Demo |
| 8.4 | Interpreting DESeq2 output | Lecture |
| **Lab 8A** | Full DESeq2 analysis on airway dataset | Hands-on |
| **Lab 8B** | Multi-factor designs in DESeq2 | Hands-on |

**Key Concepts:** GLM, dispersion estimation, results extraction

---

#### Module 9: Alternative Bulk DE Tools (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 9.1 | edgeR: empirical Bayes approach | Lecture |
| 9.2 | limma-voom: linear models | Lecture |
| 9.3 | When to use which tool | Discussion |
| **Lab 9** | Compare DESeq2, edgeR, and limma results | Hands-on |

**Key Concepts:** Tool selection, result concordance

---

### Week 4: Single-Cell DE & Interpretation

#### Module 10: Single-Cell DE Challenges (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 10.1 | Why scRNA-seq DE is different | Lecture |
| 10.2 | Dropout and sparsity | Lecture |
| 10.3 | Cell-level vs sample-level inference | Lecture |
| **Lab 10** | Explore scRNA-seq count distributions | Hands-on |

**Key Concepts:** Zero-inflation, pseudoreplication, cell dependence

---

#### Module 11: Pseudobulk Analysis (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 11.1 | Aggregating cells to samples | Lecture |
| 11.2 | Applying bulk methods to pseudobulk | Lecture |
| 11.3 | When pseudobulk is appropriate | Discussion |
| **Lab 11** | Pseudobulk DE with DESeq2 on scRNA-seq data | Hands-on |

**Key Concepts:** Aggregation, sample-level replication, proper inference

---

#### Module 12: Cell-Level DE Methods (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 12.1 | Wilcoxon tests in Seurat/Scanpy | Lecture |
| 12.2 | MAST: hurdle models | Lecture |
| 12.3 | Limitations and appropriate use | Discussion |
| **Lab 12** | Compare Wilcoxon vs pseudobulk results | Hands-on |

**Key Concepts:** Nonparametric tests, hurdle models, tool selection

---

#### Module 13: Interpreting & Visualizing DE Results (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 13.1 | Volcano plots | Lecture |
| 13.2 | MA plots | Lecture |
| 13.3 | Heatmaps of DEGs | Lecture |
| 13.4 | Gene set enrichment: GO, KEGG, Reactome | Lecture |
| **Lab 13** | Create publication-quality DE visualizations | Hands-on |

**Key Concepts:** Visualization best practices, pathway analysis

---

## Assessments

### Quizzes (after each module)
- 5-10 questions covering key concepts
- Immediate feedback

### Practical Assignments

| Assignment | Module | Description |
|------------|--------|-------------|
| **A1** | 3 | Design matrix construction for 3 scenarios |
| **A2** | 6 | Multiple testing simulation and analysis |
| **A3** | 8 | Complete DESeq2 analysis with report |
| **A4** | 11-12 | Pseudobulk vs cell-level DE comparison |

### Final Project
**End-to-end DE analysis on a published dataset**
- Choose bulk OR single-cell dataset
- Perform complete DE workflow
- Create visualization report
- Write biological interpretation (1 page)

---

## Software Setup

### R Environment (Primary)

```r
# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Core packages
BiocManager::install(c(
    "DESeq2",
    "edgeR", 
    "limma",
    "airway",        # Example dataset
    "org.Hs.eg.db",  # Annotation
    "clusterProfiler"  # Enrichment
))

# Visualization
install.packages(c("ggplot2", "pheatmap", "EnhancedVolcano"))
```

### Python Alternative

```bash
pip install scanpy anndata pydeseq2 gseapy matplotlib seaborn
```

---

## Demo Datasets

| Dataset | Type | Source | Use |
|---------|------|--------|-----|
| airway | Bulk RNA-seq | Bioconductor | Labs 4-9 |
| pasilla | Bulk RNA-seq | Bioconductor | Assignment 3 |
| PBMC 3k | scRNA-seq | 10x/Scanpy | Labs 10-12 |

---

## Resources

### Primary References
- DESeq2 vignette (Bioconductor)
- edgeR user guide
- Love MI et al. (2014) DESeq2 paper
- Soneson & Robinson (2018) - scRNA-seq DE benchmarking

### Tools Documentation
| Tool | URL |
|------|-----|
| DESeq2 | https://bioconductor.org/packages/DESeq2 |
| edgeR | https://bioconductor.org/packages/edgeR |
| limma | https://bioconductor.org/packages/limma |
| Scanpy | https://scanpy.readthedocs.io |
| clusterProfiler | https://bioconductor.org/packages/clusterProfiler |

---

## Common Mistakes to Avoid

| Mistake | Consequence | Prevention |
|---------|-------------|------------|
| No biological replicates | Cannot estimate variance | Design with ≥3 replicates |
| Ignoring batch effects | False positives/negatives | Include batch in model |
| Cell-level DE with samples | Pseudoreplication | Use pseudobulk |
| Treating p-values as truth | Overconfident conclusions | Consider effect size + biology |
| Overinterpreting small log2FC | Biologically irrelevant DEGs | Use log2FC cutoffs |
| No multiple testing correction | Massive false positives | Always use FDR |

---

## What Bad DE Analysis Looks Like (and How to Catch It)

### Symptoms (Red Flags)

- **Design**: no replicates; batch perfectly confounded with condition; unclear contrasts.
- **Normalization**: samples separate purely by library size; MA plot shows global shifts; size factors look extreme.
- **Results**: thousands of DEGs with tiny effect sizes; volcano plot dominated by very low-count genes; inconsistent direction across replicates.
- **Interpretation**: reporting “significant” genes without effect size; claiming causality from association; ignoring known biology/controls.

### Common Failure Modes → What It Looks Like → Fix

- **Pseudoreplication (treating cells as replicates)**  
  - **Looks like**: huge significance everywhere; unstable results across donors.  
  - **Fix**: use **pseudobulk** (Module 11) or sample-aware models; report sample size and donor structure.

- **Confounded design (batch = condition)**  
  - **Looks like**: perfect separation but no way to attribute to biology.  
  - **Fix**: redesign; if impossible, be explicit that batch cannot be separated; avoid strong biological claims.

- **Low-count driven DE**  
  - **Looks like**: top DEGs are barely expressed; log2FC unstable; many zeros.  
  - **Fix**: filtering, independent filtering, shrinkage, and minimum expression criteria.

- **Multiple testing ignored**  
  - **Looks like**: “significant” genes explode with p<0.05 threshold.  
  - **Fix**: always report **FDR-adjusted** p-values (BH), specify threshold (e.g., FDR < 0.05).

- **Effect size ignored**  
  - **Looks like**: statistically significant but biologically trivial results.  
  - **Fix**: require meaningful log2FC (context-specific) + baseline expression checks.

### Minimum “Sanity Checks” Before You Report

- **Replicates**: you can name the biological replicate unit and show \(n\) per group.
- **QC**: outlier samples identified and justified (not silently removed).
- **Plots**: MA + volcano + sample clustering, and they tell a consistent story.
- **Reproducibility**: script/Rmd re-runs end-to-end and generates the same tables/figures.

---

## Learning Path Progression

```
BEGINNER (Weeks 1-2)
│
├── Understand count distributions
├── Run t-tests on expression data
├── Apply FDR correction
└── Basic DESeq2 workflow
        │
        ▼
INTERMEDIATE (Weeks 3-4)
│
├── Multi-factor DE designs
├── Pseudobulk scRNA-seq
├── Compare multiple DE tools
└── Pathway enrichment analysis
        │
        ▼
ADVANCED (Post-course)
│
├── Mixed-effects models
├── Bayesian DE methods
├── Time-series DE
└── Multi-condition contrasts
```

---

## Connection to Other Courses

| Prerequisite for | Why |
|------------------|-----|
| Cell-Cell Communication | Requires identifying DEGs as ligands/receptors |
| Trajectory Analysis | DEGs along pseudotime |
| Machine Learning | Feature selection from DEGs |
| Integration Methods | Batch-corrected DE |

| Prerequisite from | Why |
|-------------------|-----|
| scRNA-seq Processing | Need clean count matrix |

---

## Suggested Schedule

| Day | Activity | Time |
|-----|----------|------|
| Mon/Wed/Fri | Lectures + Reading | 1 hr |
| Tue/Thu | Labs | 1.5 hrs |
| Weekend | Assignments | 2 hrs |

*Total: ~6 hours/week for 4 weeks*

---

## Certificate Criteria

To complete this course:
- [ ] Pass all module quizzes (≥70%)
- [ ] Submit all 4 practical assignments
- [ ] Complete final project with proper interpretation
- [ ] Document reproducible analysis

---

*This course bridges statistical theory and practical application, preparing you for rigorous transcriptomics analysis.*

