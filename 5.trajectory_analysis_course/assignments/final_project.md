# Final Project: Complete Trajectory & Fate Mapping Analysis

**Duration:** 1 week
**Points:** 250

> **Updated February 2026** — Now includes CellRank 2 multiview fate mapping requirements.

---

## Overview

Perform a comprehensive trajectory analysis and fate mapping on a differentiation dataset, integrating classical trajectory inference methods with CellRank 2's multiview framework.

---

## Dataset Options

### Option A: Hematopoiesis
- Bone marrow differentiation
- Multiple lineages (erythroid, myeloid, lymphoid)
- Rich in known markers for validation

### Option B: Organoid Development
- iPSC differentiation
- Single lineage with stages
- Well-characterized transcriptional programs

### Option C: Immune Activation
- T cell activation time course
- Cyclic/activation trajectory
- Known activation markers

### Option D: Custom Dataset
- Must have continuous differentiation process
- Must support at least 2 CellRank kernels
- Describe the system and justify dataset choice (include in report)

---

## Requirements

### 1. Data Preparation (25 points)
- QC and preprocessing
- Cell type annotation
- Verify trajectory-like structure exists (PAGA or similar)
- Describe the expected biological trajectory

### 2. Pseudotime Inference (35 points)
- Apply **2+ methods** (e.g., DPT + Monocle 3 or Slingshot)
- Compare pseudotime orderings (Spearman correlation)
- Identify root and terminal states
- Justify best ordering

### 3. RNA Velocity (30 points)
- Compute velocity (stochastic or dynamical model)
- Validate direction matches pseudotime
- Interpret velocity patterns
- Assess velocity quality (confidence scores, model fit)

### 4. CellRank 2 Multiview Fate Mapping (50 points)

> This is the new core component.

- **Run at least 2 different kernels** (15 pts)
  - VelocityKernel, CytoTRACEKernel, PseudotimeKernel, and/or RealTimeKernel
  - For each: GPCCA macrostates, terminal states, fate probabilities
- **Compare kernels** (10 pts)
  - Terminal state agreement
  - Fate probability correlation
  - Disagreement analysis
- **Combine kernels** (10 pts)
  - Weighted combination with justification
  - Compare combined vs individual results
- **Driver genes** (10 pts)
  - Identify driver genes per lineage from combined kernel
  - Validate against known markers
- **Gene trends** (5 pts)
  - Visualize gene expression trends along lineage fate probabilities

### 5. Trajectory DE (30 points)
- Identify trajectory-associated genes (pseudotime-based)
- Categorize expression patterns (early, transitional, late)
- Compare with CellRank driver genes
- Visualize gene dynamics

### 6. Biological Interpretation & Report (80 points)
- Enrichment analysis for driver/trajectory genes (15 pts)
- Connect findings to known biology (15 pts)
- **Comprehensive report (800+ words)** covering: (30 pts)
  - Introduction to the biological system
  - Methods: what you did and why
  - Results: key findings with figures
  - Discussion: biological interpretation, kernel agreement, limitations
- **Figures** — publication-quality visualizations (10 pts)
- **Validation strategy** — how would you experimentally validate? (10 pts)

---

## Deliverable Structure

```
final_project/
├── data/
├── notebooks/
│   ├── 01_preprocessing.ipynb
│   ├── 02_pseudotime.ipynb
│   ├── 03_velocity.ipynb
│   ├── 04_cellrank_individual_kernels.ipynb
│   ├── 05_cellrank_combination.ipynb
│   ├── 06_trajectory_de.ipynb
│   └── 07_gene_trends_visualization.ipynb
├── figures/
│   ├── umap_celltypes.png
│   ├── pseudotime_comparison.png
│   ├── velocity_stream.png
│   ├── cellrank_fate_per_kernel.png
│   ├── kernel_comparison.png
│   ├── combined_fate_probabilities.png
│   ├── driver_genes.png
│   └── gene_trends.png
├── results/
│   ├── pseudotime.csv
│   ├── velocity_genes.csv
│   ├── trajectory_de_genes.csv
│   ├── fate_probabilities.csv
│   ├── terminal_states.csv
│   └── driver_genes.csv
└── report.md
```

---

## Rubric

| Component | Points | Weight |
|-----------|--------|--------|
| Data preparation | 25 | 10% |
| Pseudotime inference | 35 | 14% |
| RNA velocity | 30 | 12% |
| **CellRank 2 multiview** | **50** | **20%** |
| Trajectory DE | 30 | 12% |
| Report & interpretation | 80 | 32% |
| **Total** | **250** | **100%** |

---

## Grading Scale

| Grade | Points |
|-------|--------|
| A | 225-250 |
| B | 200-224 |
| C | 175-199 |
| D | 150-174 |
| F | <150 |

---

## Tips

1. **Start with CellRank 2** — it integrates information from velocity and pseudotime, so having the CellRank analysis helps interpret earlier results
2. **Read the guide** — `resources/cellrank2_multiview_guide.md` has the complete API reference
3. **Use the CellRank Protocol** as reference — [theislab.github.io/cellrank_protocol](https://theislab.github.io/cellrank_protocol/)
4. **Validate everything** — check terminal states against known markers; check driver genes against literature
5. **Kernel disagreement is a feature** — report and interpret it, don't hide it
6. **Quality over quantity** — a thorough 2-kernel analysis beats a superficial 4-kernel analysis
7. **Cite CellRank properly** — see `resources/cellrank2_multiview_guide.md` for citation format

---

## Key References

- Weiler & Theis (2026) "CellRank protocol" *Nat. Protocols* [doi:10.1038/s41596-025-01314-w](https://doi.org/10.1038/s41596-025-01314-w)
- Weiler et al. (2024) "CellRank 2" *Nat. Methods* 21:1196-1205
- Lange et al. (2022) "CellRank for directed single-cell fate mapping" *Nat. Methods* 19:159-170
