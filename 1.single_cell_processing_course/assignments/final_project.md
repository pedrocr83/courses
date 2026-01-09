# Final Project: End-to-End scRNA-seq Processing

**Duration:** 1 week  
**Points:** 200  
**Weight:** 40% of course grade

---

## Overview

Process a scRNA-seq dataset from raw FASTQ files to a quality-controlled count matrix. Document every step for reproducibility.

---

## Dataset Options

Choose ONE of the following:

### Option A: 5k PBMCs (Recommended)
- Source: 10x Genomics
- ~5,000 cells
- Download: See setup.md

### Option B: Your Own Dataset
- Must be 10x Genomics v3 chemistry
- Provide FASTQ files
- Discuss with instructor before starting

---

## Requirements

### 1. Reference Preparation (20 points)
- Download appropriate reference genome
- Document version and source
- Build index if needed

### 2. Alignment & Quantification (40 points)
- Run Cell Ranger (or alternative pipeline)
- Document all parameters
- Record runtime and resources
- Verify successful completion

### 3. Quality Control (50 points)
- Calculate all standard QC metrics
- Generate comprehensive visualizations
- Select and justify filtering thresholds
- Run doublet detection
- Apply filters

### 4. Documentation (50 points)
- Complete the processing log template
- Complete the QC report template
- Ensure all steps are reproducible
- Include all commands and parameters

### 5. Final Deliverables (40 points)
- Filtered H5AD file
- Raw counts (unfiltered) backup
- QC summary statistics CSV
- All visualizations as PNG/PDF

---

## Submission Checklist

```
final_project/
├── data/
│   └── (symlink to FASTQ location)
├── reference/
│   └── (reference files or symlinks)
├── outputs/
│   ├── cellranger/
│   │   └── (Cell Ranger outputs)
│   ├── filtered_matrix.h5ad
│   └── qc_metrics.csv
├── figures/
│   ├── violin_plots.png
│   ├── scatter_plots.png
│   ├── barcode_rank.png
│   └── doublet_scores.png
├── docs/
│   ├── processing_log.md
│   └── qc_report.md
├── code/
│   └── qc_analysis.ipynb
└── README.md
```

---

## Evaluation Rubric

| Component | Points | Criteria |
|-----------|--------|----------|
| Reference setup | 20 | Correct version, documented source |
| Pipeline execution | 40 | Successful run, proper parameters |
| QC metrics | 25 | All metrics calculated correctly |
| Visualizations | 25 | Clear, complete, well-labeled |
| Filtering decisions | 25 | Justified, appropriate thresholds |
| Doublet handling | 15 | Detection run, results interpreted |
| Processing log | 25 | Complete, reproducible |
| QC report | 25 | Thorough, professional |

---

## Grading Scale

| Grade | Points | Description |
|-------|--------|-------------|
| A | 180-200 | Exceptional work, publication-ready |
| B | 160-179 | Strong work, minor improvements needed |
| C | 140-159 | Acceptable, some gaps in documentation |
| D | 120-139 | Incomplete or significant issues |
| F | <120 | Major components missing |

---

## Tips for Success

1. **Start early** - Cell Ranger takes hours to run
2. **Document as you go** - Don't try to remember parameters later
3. **Check intermediate outputs** - Don't wait until the end to verify
4. **Use templates** - They're designed to ensure completeness
5. **Ask questions** - Better to clarify early than redo work

---

## Academic Integrity

- All work must be your own
- You may discuss approaches with classmates
- Code must be written independently
- Cite any external resources used

