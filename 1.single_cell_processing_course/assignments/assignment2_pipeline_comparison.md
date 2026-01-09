# Assignment 2: Pipeline Comparison

**Module:** 5 - Alignment & Quantification Pipelines  
**Due:** End of Week 2  
**Points:** 100

---

## Overview

Compare Cell Ranger and kallisto|bustools pipelines by processing the same dataset with both and analyzing the differences.

---

## Tasks

### Part A: Run Both Pipelines (40 points)

1. Process the 1k PBMC dataset with Cell Ranger
2. Process the same dataset with kallisto|bustools
3. Document:
   - Commands used
   - Runtime for each
   - Resource usage (memory, CPU)

### Part B: Compare Outputs (40 points)

1. Load both count matrices into Python
2. Compare:
   - Number of cells detected
   - Number of genes detected
   - Total UMI counts
   - Median genes per cell
   - Correlation of gene expression between pipelines

3. Create visualizations:
   - Scatter plot of UMIs per cell (CellRanger vs kallisto)
   - Venn diagram of detected cells

### Part C: Discussion (20 points)

Answer the following questions (200-300 words total):

1. Which pipeline detected more cells? Why might this differ?
2. How correlated are the expression values between pipelines?
3. When would you choose one pipeline over the other?
4. What factors should influence pipeline selection in a research setting?

---

## Deliverables

1. **Report (PDF or Markdown):** 2-3 pages with visualizations
2. **Code:** Jupyter notebook with full analysis
3. **Commands file:** Text file with exact commands used for each pipeline

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Successful Cell Ranger run | 15 |
| Successful kallisto run | 15 |
| Runtime/resource documentation | 10 |
| Quantitative comparison | 20 |
| Visualizations | 15 |
| Discussion quality | 20 |
| Code quality | 5 |

---

## Tips

- Start the Cell Ranger run early—it takes time
- Use the same reference for both pipelines
- Document any errors encountered and how you resolved them

