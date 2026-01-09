# Assignment 3: Trajectory Differential Expression

**Modules:** 7-8 - Trajectory DE  
**Due:** End of Week 3  
**Points:** 100

---

## Overview

Identify genes that change along pseudotime trajectory.

---

## Tasks

### Part A: Prepare Data (20 points)

1. Use data with pseudotime from Assignment 1
2. Verify pseudotime is computed
3. Select genes for testing (expressed genes)

---

### Part B: Trajectory DE (40 points)

Using tradeSeq (R) or equivalent:

1. Fit GAM models along pseudotime
2. Test for association with pseudotime
3. Identify significantly associated genes
4. Categorize by expression pattern (increasing, decreasing, transient)

---

### Part C: Visualization (25 points)

1. Create smoothed expression plots for top genes
2. Create ordered heatmap along pseudotime
3. Identify gene modules with similar patterns

---

### Part D: Biological Interpretation (15 points)

1. Run GO enrichment on trajectory-associated genes
2. Identify key pathways activated during differentiation
3. Write 150 words on biological insights

---

## Deliverables

1. Notebook with trajectory DE
2. Figures: Gene trends, heatmap
3. Gene lists by pattern
4. Enrichment results
5. Written interpretation

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Data preparation | 20 |
| Trajectory DE | 40 |
| Visualization | 25 |
| Interpretation | 15 |

