# Final Project: Complete CCC Analysis in Disease Context

**Duration:** 1 week  
**Points:** 200  
**Weight:** 40% of course grade

---

## Overview

Conduct a complete cell-cell communication analysis on a disease-relevant dataset, from data preparation to biological interpretation.

---

## Dataset Options

### Option A: Cancer Microenvironment

- Tumor scRNA-seq with immune infiltration
- Suggested: TISCH database, GEO cancer datasets
- Focus: Tumor-immune crosstalk

### Option B: Inflammation/Autoimmune

- Inflammatory disease (IBD, RA, COVID)
- Suggested: GEO COVID PBMC studies
- Focus: Cytokine signaling

### Option C: Development

- Developmental atlas (organoids, embryo)
- Suggested: Human Cell Atlas
- Focus: Morphogen signaling (WNT, TGF-β, Notch)

### Option D: Your Own Data (Instructor approval)

---

## Requirements

### 1. Data Acquisition & Preparation (25 points)

- Load and QC the dataset
- Annotate cell types (use existing or perform)
- Document cell type composition
- Justify any filtering decisions

**Deliverables:**
- UMAP with annotations
- Cell type frequency table
- QC summary

---

### 2. CCC Analysis with 2+ Methods (50 points)

Run at least TWO of:
- CellPhoneDB
- CellChat
- NicheNet
- LIANA

For each method:
- Document parameters
- Report significant interactions
- Create standard visualizations

**Deliverables:**
- Interaction tables
- Dot plots / circle plots
- Method-specific outputs

---

### 3. Method Comparison (25 points)

- Quantify overlap between methods
- Identify consensus interactions
- Discuss disagreements

**Deliverables:**
- Venn diagram
- Comparison table
- Brief discussion

---

### 4. Network Analysis (25 points)

- Build CCC network
- Calculate centrality metrics
- Identify communication hubs
- Visualize network

**Deliverables:**
- Network figure
- Node metrics table

---

### 5. Biological Interpretation (50 points)

Write 600-800 words covering:

1. **Key Findings**
   - Most important interactions discovered
   - Hub cell types and their roles

2. **Disease Context**
   - How do findings relate to the disease?
   - Comparison to literature

3. **Novel Insights**
   - Unexpected interactions
   - Potential therapeutic targets

4. **Limitations**
   - Missing cell types?
   - Expression vs protein?
   - Spatial context missing?

5. **Validation Suggestions**
   - Experimental follow-up
   - Additional computational analyses

---

### 6. Reproducibility (25 points)

- Well-documented code
- Clear file organization
- Session info / environment
- Data availability

---

## Deliverable Structure

```
final_project/
├── data/
│   └── (processed data or links)
├── scripts/
│   ├── 01_preprocessing.R
│   ├── 02_cellchat.R
│   ├── 03_cellphonedb.py
│   ├── 04_comparison.R
│   └── 05_network.R
├── results/
│   ├── cellchat/
│   ├── cellphonedb/
│   └── network/
├── figures/
│   ├── umap_celltypes.pdf
│   ├── ccc_dotplot.pdf
│   ├── ccc_circleplot.pdf
│   ├── comparison_venn.pdf
│   └── network.pdf
├── report.pdf
└── README.md
```

---

## Report Structure

1. **Abstract** (100 words)
2. **Introduction** (200 words)
   - Disease/biological context
   - Why CCC is relevant
3. **Methods** (300 words)
   - Data source
   - Tools and parameters
4. **Results** (400 words)
   - Key interactions
   - Network structure
   - Method comparison
5. **Discussion** (400 words)
   - Biological interpretation
   - Limitations
   - Future directions
6. **Figures and Tables**
7. **References**

Total: ~1500 words + figures

---

## Evaluation Rubric

| Component | Points | Criteria |
|-----------|--------|----------|
| Data preparation | 25 | Proper QC, annotation |
| CCC analysis | 50 | 2+ methods correctly run |
| Method comparison | 25 | Systematic, insightful |
| Network analysis | 25 | Correct metrics, visualization |
| Interpretation | 50 | Biological depth, literature context |
| Reproducibility | 25 | Clear code, documented |

---

## Grading Scale

| Grade | Points | Description |
|-------|--------|-------------|
| A | 180-200 | Publication-quality analysis |
| B | 160-179 | Strong work, minor gaps |
| C | 140-159 | Adequate, some issues |
| D | 120-139 | Incomplete |
| F | <120 | Major components missing |

---

## Tips for Success

1. **Choose a focused question** - Don't try to analyze everything
2. **Start with the biology** - What interactions are you looking for?
3. **Validate computationally** - Check marker genes for your cell types
4. **Connect to literature** - Are your findings supported?
5. **Be honest about limitations** - No analysis is perfect

---

## Academic Integrity

- All work must be your own
- Cite all tools and methods
- Data must be publicly available or approved
- Reference any tutorials consulted

