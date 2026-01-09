# QC Report: [Dataset Name]

**Date:** [DATE]  
**Analyst:** [NAME]  
**Pipeline:** [Cell Ranger / kallisto / STARsolo]  
**Version:** [X.X.X]

---

## 1. Dataset Overview

| Metric | Value |
|--------|-------|
| Sample ID | |
| Species | |
| Tissue | |
| Chemistry | |
| Total cells (raw) | |
| Total cells (filtered) | |

---

## 2. Alignment Summary

| Metric | Value |
|--------|-------|
| Total reads | |
| Mapped reads (%) | |
| Reads in cells (%) | |
| Median genes per cell | |
| Median UMIs per cell | |

---

## 3. QC Thresholds Applied

| Filter | Threshold | Cells Removed |
|--------|-----------|---------------|
| Min genes | | |
| Max genes | | |
| Max mito % | | |
| Doublet score | | |

**Rationale for thresholds:**

[Explain why these thresholds were chosen]

---

## 4. QC Visualizations

### 4.1 Violin Plots
[Insert violin plots: n_genes, n_counts, pct_mito]

### 4.2 Scatter Plots
[Insert scatter: genes vs counts, colored by mito%]

### 4.3 Barcode Rank Plot
[Insert knee plot if applicable]

---

## 5. Doublet Detection

| Metric | Value |
|--------|-------|
| Method used | [Scrublet / DoubletFinder] |
| Predicted doublets | |
| Doublet rate (%) | |

---

## 6. Final Filtered Matrix

| Metric | Value |
|--------|-------|
| Cells retained | |
| Genes detected | |
| Sparsity (%) | |

---

## 7. Output Files

| File | Path | Description |
|------|------|-------------|
| Raw matrix | | Unfiltered counts |
| Filtered matrix | | After QC |
| QC metrics | | Per-cell metrics CSV |

---

## 8. Notes & Observations

[Any unusual patterns, concerns, or recommendations for downstream analysis]

---

## 9. Reproducibility

**Software versions:**
```
cellranger: 
scanpy: 
scrublet: 
python: 
```

**Reference genome:**
- Source: 
- Version: 

**Commands used:**
```bash
# Pipeline command

# Filtering script

```

