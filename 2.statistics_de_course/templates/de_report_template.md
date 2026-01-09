# Differential Expression Report: [Dataset Name]

**Date:** [DATE]  
**Analyst:** [NAME]  
**Method:** [DESeq2 / edgeR / limma / Other]  
**Version:** [X.X.X]

---

## 1. Study Overview

| Field | Value |
|-------|-------|
| Dataset | |
| Species | |
| Comparison | [Condition A] vs [Condition B] |
| Samples per group | |
| Total genes tested | |

---

## 2. Experimental Design

### Sample Information

| Sample | Condition | Batch | Other covariates |
|--------|-----------|-------|------------------|
| | | | |

### Design Formula

```
~ [design formula used]
```

**Rationale:**

[Explain why this design was chosen]

---

## 3. Quality Control

### Pre-filtering
- Genes removed (low counts): 
- Genes retained: 

### Sample QC
- PCA plot: [Reference figure]
- Sample correlation: [Reference figure]
- Outliers identified: 

---

## 4. Normalization

| Method | Details |
|--------|---------|
| Normalization | [median-of-ratios / TMM / other] |
| Size factors range | min - max |

---

## 5. Differential Expression Results

### Summary Statistics

| Metric | Value |
|--------|-------|
| Total DEGs (FDR < 0.05) | |
| Upregulated | |
| Downregulated | |
| DEGs with \|log2FC\| > 1 | |

### Thresholds Used

| Parameter | Value | Justification |
|-----------|-------|---------------|
| FDR cutoff | | |
| log2FC cutoff | | |
| baseMean cutoff | | |

---

## 6. Visualizations

### 6.1 Volcano Plot
[Insert or reference volcano plot]

### 6.2 MA Plot
[Insert or reference MA plot]

### 6.3 Heatmap of Top DEGs
[Insert or reference heatmap]

### 6.4 PCA of Samples
[Insert or reference PCA]

---

## 7. Top Differentially Expressed Genes

### Top 10 Upregulated

| Gene | log2FC | padj | baseMean |
|------|--------|------|----------|
| | | | |

### Top 10 Downregulated

| Gene | log2FC | padj | baseMean |
|------|--------|------|----------|
| | | | |

---

## 8. Pathway Enrichment

### Method Used
[GO / KEGG / Reactome / GSEA]

### Top Enriched Pathways (Upregulated genes)

| Pathway | p.adjust | Gene count |
|---------|----------|------------|
| | | |

### Top Enriched Pathways (Downregulated genes)

| Pathway | p.adjust | Gene count |
|---------|----------|------------|
| | | |

---

## 9. Biological Interpretation

[2-3 paragraphs discussing:
- What the results suggest biologically
- Key pathways and their relevance
- Comparison to published literature
- Limitations and caveats]

---

## 10. Reproducibility

### Software Versions

```
R: 
DESeq2: 
ggplot2: 
clusterProfiler: 
```

### Code Location

```
[Path to analysis scripts]
```

### Key Commands

```r
# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(...)
dds <- DESeq(dds)
res <- results(dds, ...)
```

---

## 11. Output Files

| File | Description | Location |
|------|-------------|----------|
| Full results table | All genes with statistics | |
| Significant DEGs | Filtered DEG list | |
| Normalized counts | For downstream use | |
| Figures | All visualizations | |

---

## 12. Conclusions

[Summary of key findings and recommendations for follow-up]

---

## 13. References

[Relevant papers and methods citations]

