# Final Project: Complete Trajectory Analysis

**Duration:** 1 week  
**Points:** 200

---

## Overview

Perform comprehensive trajectory analysis on a differentiation dataset.

---

## Dataset Options

### Option A: Hematopoiesis
- Bone marrow differentiation
- Multiple lineages

### Option B: Organoid Development
- iPSC differentiation
- Single lineage with stages

### Option C: Immune Activation
- T cell activation time course
- Cyclic/activation trajectory

---

## Requirements

### 1. Data Preparation (25 points)
- QC and preprocessing
- Cell type annotation
- Verify trajectory-like structure exists

### 2. Pseudotime Inference (40 points)
- Apply 2+ methods
- Compare and justify best ordering
- Identify root and terminal states

### 3. RNA Velocity (30 points)
- Compute velocity
- Validate direction matches pseudotime
- Interpret velocity patterns

### 4. Trajectory DE (35 points)
- Identify trajectory-associated genes
- Categorize expression patterns
- Visualize gene dynamics

### 5. CellRank Analysis (30 points)
- Compute fate probabilities
- Identify terminal states
- Find driver genes

### 6. Biological Interpretation (40 points)
- Enrichment analysis
- Connect to known biology
- Write comprehensive report (600 words)

---

## Deliverable Structure

```
final_project/
├── data/
├── notebooks/
│   ├── 01_preprocessing.ipynb
│   ├── 02_pseudotime.ipynb
│   ├── 03_velocity.ipynb
│   ├── 04_trajectory_de.ipynb
│   └── 05_cellrank.ipynb
├── figures/
├── results/
│   ├── pseudotime.csv
│   ├── velocity_genes.csv
│   ├── trajectory_de_genes.csv
│   └── fate_probabilities.csv
└── report.md
```

---

## Rubric

| Component | Points |
|-----------|--------|
| Data preparation | 25 |
| Pseudotime | 40 |
| Velocity | 30 |
| Trajectory DE | 35 |
| CellRank | 30 |
| Report | 40 |

---

## Grading Scale

| Grade | Points |
|-------|--------|
| A | 180-200 |
| B | 160-179 |
| C | 140-159 |
| D | 120-139 |
| F | <120 |

