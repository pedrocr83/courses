# Final Project: End-to-End DE Analysis

**Duration:** 1 week  
**Points:** 200  
**Weight:** 40% of course grade

---

## Overview

Conduct a complete differential expression analysis on a real-world dataset, from raw counts to biological interpretation.

---

## Dataset Options

### Option A: Bulk RNA-seq (Recommended for beginners)

Choose from GEO/SRA:
- **GSE93819** - Breast cancer treatment response
- **GSE60450** - Mouse diet study
- **GSE52778** - Inflammatory bowel disease

Or use:
- recount3 package for easy access

### Option B: Single-Cell with Pseudobulk

- Download from cellxgene: https://cellxgene.cziscience.com/
- Must have sample-level metadata (multiple donors)
- Perform pseudobulk analysis

### Option C: Your Own Data (Instructor approval required)

---

## Requirements

### 1. Data Acquisition & Description (20 points)

- Download and load dataset
- Document:
  - Source and accession
  - Biological question
  - Sample sizes
  - Experimental design

### 2. Quality Control (30 points)

- Assess data quality
- Filter genes appropriately
- Check for outliers
- Verify experimental design

**Deliverables:**
- PCA plot
- Sample correlation heatmap
- QC summary table

### 3. Normalization (20 points)

- Apply appropriate normalization
- Justify method choice
- Show effect of normalization

### 4. Differential Expression (40 points)

- Build appropriate design matrix
- Run DESeq2 (or equivalent)
- Apply FDR correction
- Extract significant genes

**Deliverables:**
- Full results table
- Summary statistics
- MA and volcano plots

### 5. Visualization (30 points)

Create publication-quality figures:
- Volcano plot with labeled genes
- Heatmap of top DEGs
- Expression plots for key genes
- (Optional) Pathway enrichment visualization

### 6. Biological Interpretation (40 points)

- Perform pathway/GO enrichment
- Interpret top pathways
- Connect findings to literature
- Discuss biological significance

Write 500-750 words.

### 7. Reproducibility (20 points)

- Well-documented code
- Session info
- Data availability statement
- Clear file organization

---

## Deliverable Structure

```
final_project/
├── data/
│   └── (raw data or link)
├── scripts/
│   └── analysis.Rmd (or .py)
├── results/
│   ├── de_results.csv
│   └── enrichment_results.csv
├── figures/
│   ├── pca_plot.pdf
│   ├── volcano_plot.pdf
│   ├── heatmap.pdf
│   └── enrichment_plot.pdf
├── report.pdf
└── README.md
```

---

## Report Structure

1. **Abstract** (100 words)
2. **Introduction** (200 words)
   - Biological background
   - Experimental question
3. **Methods** (300 words)
   - Data source
   - Analysis workflow
   - Tools and parameters
4. **Results** (400 words)
   - QC findings
   - DE summary
   - Key genes and pathways
5. **Discussion** (300 words)
   - Biological interpretation
   - Limitations
   - Future directions
6. **References**

Total: ~1500 words + figures

---

## Evaluation Rubric

| Component | Points | Criteria |
|-----------|--------|----------|
| Data handling | 20 | Correct loading, appropriate filtering |
| QC analysis | 30 | Thorough, well-visualized |
| DE analysis | 40 | Correct workflow, proper statistics |
| Visualization | 30 | Publication-quality, informative |
| Interpretation | 40 | Biologically meaningful, literature-informed |
| Reproducibility | 20 | Clear code, documented, reproducible |
| Report quality | 20 | Well-written, organized, complete |

---

## Grading Scale

| Grade | Points | Description |
|-------|--------|-------------|
| A | 180-200 | Exceptional, publication-ready |
| B | 160-179 | Strong analysis, minor gaps |
| C | 140-159 | Adequate, some issues |
| D | 120-139 | Incomplete or significant problems |
| F | <120 | Major components missing |

---

## Timeline

| Day | Milestone |
|-----|-----------|
| 1-2 | Data acquisition, QC |
| 3-4 | DE analysis |
| 5 | Enrichment, visualization |
| 6 | Report writing |
| 7 | Final review, submission |

---

## Tips for Success

1. **Start with the question** - What comparison are you making?
2. **Document everything** - Future you will thank present you
3. **Don't over-interpret** - Be honest about limitations
4. **Use the templates** - They ensure completeness
5. **Ask questions early** - Don't wait until day 6

---

## Academic Integrity

- All work must be your own
- Cite all software and methods
- Reference any consulted tutorials
- Data must be publicly available or approved

