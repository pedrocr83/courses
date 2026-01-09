# Statistics & Differential Expression (DE) for Transcriptomics

A structured study course for **complete beginners**, designed to build intuition first, then rigor, and finally practical mastery. This document is suitable as the backbone of a self-study course, onboarding material, or curriculum for computational biology.

---

## 1. What Is Transcriptomics?

Transcriptomics studies **RNA abundance** to understand which genes are active, where, and under what conditions.

### Common Data Types
- **Bulk RNA-seq**: average expression across many cells
- **Single-cell RNA-seq (scRNA-seq)**: expression per cell
- **Spatial transcriptomics**: expression with spatial coordinates

### Core Goal
Compare gene expression between conditions, cell types, or states to identify **differentially expressed genes (DEGs)**.

---

## 2. Why Statistics Is Central to Transcriptomics

Gene expression measurements are:
- Noisy
- Sparse (especially scRNA-seq)
- High-dimensional (10k–30k genes)

Statistics answers:
- Is a gene truly different or just noisy?
- How confident are we?
- How many false positives are acceptable?

---

## 3. Data Representation Basics

### Count Matrix
- Rows: genes
- Columns: samples or cells
- Values: read/UMI counts

### Properties of Count Data
- Non-negative integers
- Skewed distributions
- Mean–variance dependence

---

## 4. Probability & Distributions (Minimal but Essential)

### Distributions You Must Know

#### Poisson Distribution
- Models random counts
- Assumes mean = variance
- Too simple for RNA-seq

#### Negative Binomial Distribution
- Allows overdispersion
- Backbone of bulk RNA-seq DE
- Used by DESeq2 and edgeR

#### Zero-Inflated Models
- Handle excess zeros
- Relevant for scRNA-seq

---

## 5. Experimental Design (Critical Section)

Bad design cannot be fixed statistically.

### Key Concepts
- Biological replicates vs technical replicates
- Confounders (batch, donor, sex, protocol)
- Balanced designs

### Design Matrix Intuition
A table encoding which variables explain expression variation.

---

## 6. Normalization: Making Samples Comparable

### Why Normalize?
Sequencing depth and capture efficiency vary.

### Common Methods

#### Bulk RNA-seq
- Library size normalization
- TMM (edgeR)
- Median-of-ratios (DESeq2)

#### Single-cell RNA-seq
- Counts per 10k (CP10K)
- Log-normalization
- SCTransform (variance stabilization)

---

## 7. Differential Expression: Core Idea

**Differential expression** asks:
> Is the expression of gene X systematically different between groups A and B?

Statistical testing is performed **per gene**.

---

## 8. Hypothesis Testing Fundamentals

### Null Hypothesis (H₀)
No difference in expression between groups.

### Test Statistic
A number summarizing evidence against H₀.

### p-value
Probability of observing the data (or more extreme) if H₀ is true.

Important: p-values are *not* effect sizes.

---

## 9. Multiple Testing Problem

Testing 20,000 genes → many false positives.

### Corrections
- Bonferroni (too strict)
- False Discovery Rate (FDR)

### Benjamini–Hochberg (Standard)
Controls expected proportion of false positives.

---

## 10. Effect Size & Biological Meaning

Statistical significance ≠ biological relevance.

### Common Effect Sizes
- Log2 fold change (log2FC)
- Average expression difference

Always interpret DEGs using:
- Effect size
- Expression level
- Known biology

---

## 11. Bulk RNA-seq DE Methods

### DESeq2
- Negative binomial GLM
- Shrinkage of fold changes
- Gold standard for bulk RNA-seq

Docs: https://bioconductor.org/packages/release/bioc/html/DESeq2.html

### edgeR
- Empirical Bayes estimation
- Flexible experimental designs

Docs: https://bioconductor.org/packages/release/bioc/html/edgeR.html

### limma-voom
- Linear models after variance modeling
- Fast, robust

Docs: https://bioconductor.org/packages/release/bioc/html/limma.html

---

## 12. Single-Cell Differential Expression

### Why It’s Harder
- Dropout
- Cell–cell dependence
- Large sample sizes

### Common Approaches

#### Pseudobulk (Recommended)
- Aggregate cells per sample
- Apply bulk RNA-seq DE

#### Cell-Level Tests
- Wilcoxon rank-sum (Scanpy/Seurat)
- Model-based (MAST)

---

## 13. Popular scRNA-seq DE Tools

### Seurat
- Wilcoxon, logistic regression

Docs: https://satijalab.org/seurat/

### Scanpy
- Nonparametric tests

Docs: https://scanpy.readthedocs.io/

### MAST
- Hurdle model for scRNA-seq

Paper: https://pubmed.ncbi.nlm.nih.gov/26653891/

---

## 14. Interpreting DE Results

### Volcano Plots
- Effect size vs significance

### Heatmaps
- Expression patterns across samples

### Gene Set Enrichment
- GO
- KEGG
- Reactome

DE is a starting point, not the end.

---

## 15. Common Pitfalls

- Ignoring batch effects
- Overinterpreting small log2FC
- Using cell-level DE when samples exist
- Treating p-values as truth

---

## 16. From Statistics to Biology

Ask:
- Which pathways change?
- Which cell types drive differences?
- Are results reproducible?

DE supports hypotheses; it does not prove mechanisms.

---

## 17. Learning Path (Beginner → Advanced)

### Beginner
- Understand distributions
- Run DESeq2 on toy data

### Intermediate
- Design matrices
- Pseudobulk scRNA-seq

### Advanced
- Mixed models
- Bayesian DE
- Multi-condition contrasts

---

## 18. Preparing for Advanced Topics

This document prepares you for:
- Cell–cell communication inference
- Trajectory analysis
- Perturbation modeling
- ML-based expression modeling

---

## 19. Suggested Outcome

By the end, you should be able to:
- Understand RNA-seq statistics
- Perform DE analysis correctly
- Interpret results biologically
- Avoid common analytical traps

---

*This document is intentionally written for clarity over mathematical density and can be extended into notebooks, lectures, or exams.*

