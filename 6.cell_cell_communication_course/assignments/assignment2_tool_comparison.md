# Assignment 2: CellPhoneDB vs CellChat Comparison

**Modules:** 4-5 - Core CCC Tools  
**Due:** End of Week 2  
**Points:** 150

---

## Overview

Run both CellPhoneDB and CellChat on the same dataset and systematically compare their results.

---

## Dataset

Use your prepared PBMC data from Assignment 1.

---

## Tasks

### Part A: Run CellPhoneDB (35 points)

1. Run CellPhoneDB with default parameters
2. Document:
   - Runtime
   - Number of significant interactions (p < 0.05)
   - Top 10 interactions by significance
3. Create visualizations:
   - Dot plot of top interactions
   - Heatmap of interaction counts

```python
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# Run analysis
results = cpdb_statistical_analysis_method.call(
    cpdb_file_path=cpdb_db_path,
    meta_file_path='meta.txt',
    counts_file_path='counts.txt',
    counts_data='hgnc_symbol',
    output_path='results/'
)
```

---

### Part B: Run CellChat (35 points)

1. Run CellChat with default parameters
2. Document:
   - Runtime
   - Number of significant interactions
   - Top 10 interactions by probability
3. Create visualizations:
   - Circle plot of overall communication
   - Pathway-level heatmap

```r
library(CellChat)

# Create object
cellchat <- createCellChat(object = seurat_obj, group.by = "cell_type")

# Set database
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# Run analysis
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
```

---

### Part C: Systematic Comparison (50 points)

1. **Quantitative Comparison:**
   - How many interactions does each tool identify?
   - How many are in common (same L-R, same cell types)?
   - Calculate Jaccard similarity

2. **Venn Diagram:**
   - Create Venn diagram of significant interactions

3. **Top Interactions Overlap:**
   - Do the top 20 interactions agree?
   - List agreements and disagreements

4. **Pathway Comparison (CellChat only):**
   - Which pathways are most active?
   - Do CellPhoneDB results support the same pathways?

5. **Scoring Comparison:**
   - For shared interactions, is the ranking similar?
   - Scatter plot: CellPhoneDB score vs CellChat probability

---

### Part D: Critical Evaluation (30 points)

Write 300-400 words addressing:

1. **Why do the tools disagree?**
   - Different databases?
   - Different statistical methods?
   - Different scoring approaches?

2. **Which results do you trust more? Why?**

3. **When would you use each tool?**
   - CellPhoneDB strengths
   - CellChat strengths

4. **What would you report in a publication?**

---

## Deliverables

1. **Code:** Notebooks/scripts for both tools
2. **Figures:**
   - CellPhoneDB dot plot
   - CellChat circle plot
   - Venn diagram of overlap
   - Correlation scatter plot
3. **Tables:**
   - Top 20 interactions from each tool
   - Overlap summary statistics
4. **Written evaluation**

---

## Rubric

| Criterion | Points |
|-----------|--------|
| CellPhoneDB analysis complete | 35 |
| CellChat analysis complete | 35 |
| Systematic comparison thorough | 50 |
| Critical evaluation insightful | 30 |

---

## Expected Findings

- Tools will identify different numbers of interactions
- ~30-60% overlap is typical
- Different databases contribute to differences
- Neither tool is "correct" - they have different biases

