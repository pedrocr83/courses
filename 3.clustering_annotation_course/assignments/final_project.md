# Final Project: Annotate a Novel Dataset

**Duration:** 1 week  
**Points:** 200

---

## Overview

Perform complete clustering and annotation analysis on a dataset you haven't seen before.

---

## Dataset Options

### Option A: PBMC 10k
- More cells, more complexity than 3k
- 10x Genomics website

### Option B: Tissue Dataset
- Choose from cellxgene (non-PBMC tissue)
- Must have >5,000 cells

### Option C: Instructor-Provided
- Mystery dataset for testing

---

## Requirements

### 1. Preprocessing (20 points)
- QC filtering
- Normalization
- HVG selection
- PCA

### 2. Dimensionality Reduction (30 points)
- PCA analysis with elbow plot
- UMAP with parameter exploration
- Justified parameter choices

### 3. Clustering (40 points)
- Test multiple resolutions
- Stability analysis
- Justified resolution choice

### 4. Marker Identification (30 points)
- DE for markers
- Visualizations (dot plot, heatmap)

### 5. Manual Annotation (40 points)
- Use appropriate marker databases
- Document evidence for each annotation
- Handle ambiguous clusters

### 6. Automated Annotation (20 points)
- Run 1+ automated method
- Compare to manual

### 7. Final Deliverable (20 points)
- Publication-quality UMAP
- Complete annotation table
- Brief report (500 words)

---

## Deliverable Structure

```
final_project/
├── data/
├── notebooks/
│   └── analysis.ipynb
├── figures/
│   ├── qc_plots.png
│   ├── pca_elbow.png
│   ├── umap_clusters.png
│   ├── marker_dotplot.png
│   └── final_annotated_umap.png
├── results/
│   ├── markers.csv
│   └── annotations.csv
└── report.md
```

---

## Rubric

| Component | Points |
|-----------|--------|
| Preprocessing | 20 |
| Dimensionality reduction | 30 |
| Clustering | 40 |
| Markers | 30 |
| Manual annotation | 40 |
| Automated annotation | 20 |
| Final deliverables | 20 |

---

## Grading Scale

| Grade | Points |
|-------|--------|
| A | 180-200 |
| B | 160-179 |
| C | 140-159 |
| D | 120-139 |
| F | <120 |

