# Final Project: Multi-Batch Atlas Construction

**Duration:** 1 week  
**Points:** 200

---

## Overview

Build an integrated atlas from multiple scRNA-seq batches and perform downstream analysis.

---

## Dataset Options

### Option A: Multi-Donor PBMC
- 6-8 donors
- Identify batch vs biological variation

### Option B: Cross-Technology
- Same tissue, different platforms
- Available from benchmarking studies

### Option C: Multi-Condition
- Control vs disease
- Integrate while preserving condition differences

---

## Requirements

### 1. Data Preparation (25 points)
- QC each batch separately
- Document batch properties
- Identify potential confounding

### 2. Batch Effect Diagnosis (25 points)
- Visual assessment
- Quantitative metrics pre-integration

### 3. Integration (50 points)
- Apply 2+ integration methods
- Document parameters
- Justify method choice

### 4. Evaluation (40 points)
- Full metric evaluation
- Method comparison
- Quality assessment

### 5. Downstream Analysis (40 points)

On your integrated data, perform:
- Clustering
- Cell type annotation
- Condition comparison (if applicable)

### 6. Report (20 points)
- 500-word summary
- All figures publication-quality
- Methods documentation

---

## Deliverable Structure

```
final_project/
├── data/
├── notebooks/
│   ├── 01_qc_per_batch.ipynb
│   ├── 02_integration.ipynb
│   └── 03_downstream.ipynb
├── figures/
├── results/
│   ├── integrated_adata.h5ad
│   └── metrics.csv
└── report.md
```

---

## Rubric

| Component | Points |
|-----------|--------|
| Data preparation | 25 |
| Diagnosis | 25 |
| Integration | 50 |
| Evaluation | 40 |
| Downstream | 40 |
| Report | 20 |

