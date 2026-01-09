# Assignment 3: QC Report with Filtering Decisions

**Modules:** 8-9 - Quality Control & Filtering  
**Due:** End of Week 3  
**Points:** 100

---

## Overview

Create a comprehensive QC report for a scRNA-seq dataset, document your filtering decisions, and justify your thresholds.

---

## Tasks

### Part A: QC Metrics Visualization (30 points)

Using the 1k PBMC processed data:

1. Generate violin plots for:
   - Genes per cell
   - UMIs per cell
   - Mitochondrial percentage

2. Generate scatter plots for:
   - UMIs vs genes (colored by mito%)
   - UMIs vs mito%

3. Generate a barcode rank plot (if using raw data)

### Part B: Threshold Selection (30 points)

1. Determine appropriate thresholds for:
   - Minimum genes per cell
   - Maximum genes per cell
   - Maximum mitochondrial percentage

2. For each threshold, provide:
   - The value you chose
   - How you determined it (visual inspection, MAD-based, literature)
   - How many cells would be removed

3. Compare fixed thresholds vs MAD-based adaptive thresholds

### Part C: Doublet Detection (20 points)

1. Run Scrublet on your filtered data
2. Report:
   - Number of predicted doublets
   - Doublet rate
   - Whether the automatic threshold seems appropriate

3. Visualize the doublet score distribution

### Part D: Final Summary (20 points)

1. Create a summary table showing:
   - Cells before/after each filtering step
   - Genes before/after filtering
   - Final matrix dimensions

2. Export your filtered matrix to H5AD format

3. Write a brief paragraph (100-150 words) discussing:
   - Overall data quality
   - Any concerns for downstream analysis
   - Recommendations for the analysis team

---

## Deliverables

1. **QC Report:** Use the template in `templates/qc_report_template.md`
2. **Jupyter Notebook:** Complete analysis code
3. **Filtered Data:** H5AD file of your filtered matrix

---

## Rubric

| Criterion | Points |
|-----------|--------|
| QC visualizations (complete & clear) | 15 |
| Threshold documentation | 15 |
| Threshold justification | 15 |
| Doublet detection | 15 |
| Before/after comparison | 10 |
| Summary discussion | 10 |
| Proper data export | 10 |
| Code quality | 10 |

---

## Tips

- Use the QC report template—it will guide your documentation
- Don't over-filter! Justify removing each cell
- Compare your thresholds to published PBMC studies

