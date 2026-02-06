# Assignment 4: CellRank 2 Multiview Fate Mapping

**Modules 7-12** — CellRank 2 Framework
**Duration:** 1 week
**Points:** 100

---

## Overview

Apply CellRank 2's multiview framework to a differentiation dataset. You will run multiple kernels, compare their results, combine them, and draw biological conclusions about cell fate.

---

## Dataset

Use **one** of the following:

### Option A: Pancreas Endocrine Differentiation
```python
import scvelo as scv
adata = scv.datasets.pancreas()
```
- Has velocity data → VelocityKernel, CytoTRACEKernel, PseudotimeKernel

### Option B: Bone Marrow Hematopoiesis
```python
adata = scv.datasets.bonemarrow()
```
- Multi-lineage differentiation → CytoTRACEKernel, PseudotimeKernel

### Option C: Your Own Dataset
- Must have at least 2 applicable kernels
- Describe the biological system in your report

---

## Requirements

### Part 1: Individual Kernel Analysis (30 points)

Run at least **2 different CellRank kernels** on your data. For each kernel:

1. **Compute the transition matrix** (5 pts)
2. **Run GPCCA** — Schur decomposition + macrostates (5 pts)
3. **Identify terminal states** — justify your `n_states` choice with the eigenvalue gap (5 pts)
4. **Compute fate probabilities** — visualize on UMAP (5 pts)
5. **Find top 10 driver genes** per lineage (5 pts)
6. **Interpret biologically** — do terminal states match known biology? (5 pts)

### Part 2: Kernel Comparison (25 points)

Compare your kernels:

1. **Terminal state comparison** — do both kernels find the same endpoints? (5 pts)
2. **Fate probability correlation** — Spearman correlation between kernels for shared lineages (5 pts)
3. **Disagreement analysis** — identify cells/clusters with highest kernel disagreement (5 pts)
4. **Biological interpretation** — explain why kernels agree or disagree (5 pts)
5. **Visualization** — side-by-side fate probability UMAPs + correlation heatmap (5 pts)

### Part 3: Kernel Combination (25 points)

Combine your kernels:

1. **Choose weights** — justify your weight selection (5 pts)
2. **Run GPCCA on combined kernel** (5 pts)
3. **Compare combined vs individual** — are terminal states more/less stable? (5 pts)
4. **Driver genes from combined analysis** — compare with individual kernel drivers (5 pts)
5. **Weight sensitivity** — try 2+ different weight combinations; report stability (5 pts)

### Part 4: Report (20 points)

Write a **500-word report** covering:

1. **Methods** — which kernels, why, what weights, what parameters (5 pts)
2. **Results** — key findings, terminal states, fate probabilities, drivers (5 pts)
3. **Agreement/disagreement** — what you learned from comparing kernels (5 pts)
4. **Biological conclusions** — what the fate mapping tells you about the system (5 pts)

---

## Deliverables

```
assignment4/
├── notebooks/
│   ├── 01_kernel_individual.ipynb     # Part 1: each kernel separately
│   ├── 02_kernel_comparison.ipynb     # Part 2: comparison analysis
│   └── 03_kernel_combination.ipynb    # Part 3: combined analysis
├── figures/
│   ├── fate_probabilities_*.png       # UMAP plots per kernel
│   ├── kernel_correlation.png         # Correlation heatmap
│   └── combined_analysis.png          # Combined kernel results
├── results/
│   ├── terminal_states.csv            # Terminal states per kernel
│   ├── driver_genes.csv               # Driver genes per kernel/lineage
│   └── fate_probabilities.csv         # Fate probability matrix
└── report.md                          # 500-word report
```

---

## Rubric

| Component | Points | Criteria |
|-----------|--------|----------|
| Individual kernels | 30 | Correct API usage; reasonable parameters; biological interpretation |
| Kernel comparison | 25 | Quantitative comparison; clear visualization; thoughtful analysis |
| Kernel combination | 25 | Justified weights; sensitivity analysis; comparison with individual |
| Report | 20 | Clear writing; correct terminology; biological insight |
| **Total** | **100** | |

### Grading Scale

| Grade | Points |
|-------|--------|
| A | 90-100 |
| B | 80-89 |
| C | 70-79 |
| D | 60-69 |
| F | <60 |

---

## Tips

- Start by reading `resources/cellrank2_multiview_guide.md` for the API reference
- The CellRank Protocol Jupyter Book has complete example notebooks: [theislab.github.io/cellrank_protocol](https://theislab.github.io/cellrank_protocol/)
- Always plot the eigenvalue spectrum before choosing `n_states`
- If kernels disagree, that's not a failure — it's a finding to investigate
- Use the same `adata` object for all kernels (essential for combination)
- Save intermediate results — CellRank computations can be slow on large datasets
