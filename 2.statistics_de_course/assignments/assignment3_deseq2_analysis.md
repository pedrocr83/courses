# Assignment 3: Complete DESeq2 Analysis

**Module:** 8 - Bulk RNA-seq DE with DESeq2  
**Due:** End of Week 3  
**Points:** 150

---

## Overview

Perform a complete differential expression analysis on the `pasilla` dataset and write a comprehensive report.

---

## Dataset: pasilla

Drosophila melanogaster RNA-seq comparing knockdown of pasilla gene vs control.

```r
BiocManager::install("pasilla")
library(pasilla)

# Load data
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv",
                       package = "pasilla", mustWork = TRUE)
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv",
                        package = "pasilla", mustWork = TRUE)

cts <- as.matrix(read.csv(pasCts, sep="\t", row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
```

---

## Part A: Data Preparation (25 points)

1. Load and inspect the count matrix and sample metadata
2. Ensure samples are in the same order
3. Create appropriate factor levels (untreated as reference)
4. Build the DESeqDataSet with appropriate design
5. Filter low-count genes (keep genes with ≥10 counts in ≥3 samples)

**Document:**
- Starting number of genes
- Number of genes after filtering
- Sample sizes per group

---

## Part B: Quality Control (25 points)

1. Run variance-stabilizing transformation (vst)
2. Create PCA plot colored by condition
3. Create sample distance heatmap
4. Assess for outliers or batch effects

**Include:**
- PCA plot
- Distance heatmap
- Brief interpretation (2-3 sentences)

---

## Part C: Differential Expression (40 points)

1. Run DESeq2 analysis
2. Extract results for treated vs untreated
3. Apply log2FC shrinkage using `lfcShrink()`
4. Summarize results at different thresholds:
   - FDR < 0.1
   - FDR < 0.05
   - FDR < 0.01

**Report:**
- Total significant DEGs at each threshold
- Number upregulated vs downregulated
- Top 10 DEGs by adjusted p-value

---

## Part D: Visualization (30 points)

Create publication-quality figures:

1. **MA Plot** - log2FC vs mean expression
2. **Volcano Plot** - log2FC vs -log10(p-value)
3. **Heatmap** - top 50 DEGs across samples
4. **Gene expression boxplots** - 2-3 genes of interest

**Requirements:**
- Clear labels and legends
- Appropriate color schemes
- Significance thresholds marked

---

## Part E: Interpretation (30 points)

Write a brief report (300-500 words) addressing:

1. Summary of findings
2. Biological interpretation (what genes/functions are affected?)
3. Any concerns or caveats
4. Suggestions for follow-up experiments

---

## Deliverables

1. **R Markdown document** with all code
2. **Knitted HTML or PDF report**
3. **Figures** (saved as PNG/PDF)
4. **Results table** (CSV of all significant DEGs)

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Data preparation correct | 25 |
| QC analysis complete | 25 |
| DESeq2 run correctly | 25 |
| Results extraction proper | 15 |
| Visualization quality | 30 |
| Written interpretation | 20 |
| Code quality & reproducibility | 10 |

---

## Starter Code

```r
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load pasilla data
# ... your code here ...

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ condition
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, alpha = 0.05)

# Shrinkage
resLFC <- lfcShrink(dds, coef = "condition_treated_vs_untreated", type = "apeglm")
```

---

## Resources

- DESeq2 vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- pasilla dataset documentation
- EnhancedVolcano package for volcano plots

