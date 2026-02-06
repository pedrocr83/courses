# CellRank 2: Multiview Fate Mapping Reference Guide

**Quick reference for CellRank 2 kernels, estimators, and workflows.**

Based on: Weiler & Theis (2026) "CellRank: consistent and data view agnostic fate mapping for single-cell genomics." *Nature Protocols*. [doi:10.1038/s41596-025-01314-w](https://doi.org/10.1038/s41596-025-01314-w)

---

## Architecture Overview

```
DATA VIEWS                    KERNELS                     ESTIMATORS
═══════════                   ═══════                     ══════════

RNA velocity      ──→  VelocityKernel     ─┐
                                            │
Gene counts       ──→  CytoTRACEKernel    ─┤
                                            ├──→  GPCCA  ──→  Macrostates
Pseudotime        ──→  PseudotimeKernel   ─┤                  Terminal states
                                            │                  Fate probabilities
Experimental time ──→  RealTimeKernel     ─┤                  Driver genes
                                            │                  Gene trends
kNN graph         ──→  ConnectivityKernel ─┘

           Kernels can be COMBINED:
           combined = 0.8 * vk + 0.2 * ck
```

---

## Kernel Reference

### VelocityKernel

**Input:** RNA velocity vectors (from scVelo)
**Biological signal:** Local transcriptional dynamics via spliced/unspliced RNA
**When to use:** When you have high-quality velocity estimates

```python
from cellrank.kernels import VelocityKernel

vk = VelocityKernel(adata)
vk.compute_transition_matrix()
```

**Prerequisites:**
- `adata.layers['spliced']` and `adata.layers['unspliced']`
- scVelo velocity computed (`scv.tl.velocity()`)
- `adata.uns['velocity_graph']`

**Strengths:** Local, gene-level directional information
**Limitations:** Depends on velocity quality; may be noisy; assumes steady-state or dynamical model

---

### CytoTRACEKernel

**Input:** Gene expression matrix (raw or normalized)
**Biological signal:** Developmental potential via gene count complexity
**When to use:** No velocity data available; developmental/differentiation systems

```python
from cellrank.kernels import CytoTRACEKernel

ctk = CytoTRACEKernel(adata)
ctk.compute_transition_matrix()

# CytoTRACE scores stored in:
# adata.obs['ct_pseudotime']  (0 = stem, 1 = differentiated)
```

**Prerequisites:** Standard preprocessed scRNA-seq data (no special layers needed)

**Strengths:** Always available; no velocity needed; robust for development
**Limitations:** Assumes gene count decreases with differentiation; may not work for all systems (e.g., T cell activation)

---

### PseudotimeKernel

**Input:** Any pseudotime stored in `adata.obs`
**Biological signal:** Whatever the pseudotime method captures
**When to use:** When you have a trusted pseudotime ordering from any method

```python
from cellrank.kernels import PseudotimeKernel

ptk = PseudotimeKernel(adata, time_key='dpt_pseudotime')
ptk.compute_transition_matrix()
```

**Compatible pseudotime sources:**
- Diffusion pseudotime (`sc.tl.dpt()`)
- Monocle 3 pseudotime
- Slingshot pseudotime
- CytoTRACE pseudotime
- Any custom ordering in `adata.obs`

**Strengths:** Wraps any pseudotime into CellRank; flexible
**Limitations:** Quality depends entirely on input pseudotime; sensitive to root cell choice

---

### RealTimeKernel

**Input:** Experimental time points in `adata.obs`
**Biological signal:** Actual measured temporal progression
**When to use:** Time-course scRNA-seq experiments

```python
from cellrank.kernels import RealTimeKernel

rtk = RealTimeKernel(adata, time_key='experimental_time')
rtk.compute_transition_matrix(threshold='auto')
```

**Prerequisites:**
- Column in `adata.obs` with experimental time points (numeric)
- Cells collected at discrete time points

**Strengths:** Uses real time (strongest directional signal); no assumptions about splicing or gene counts
**Limitations:** Only for time-course experiments; optimal transport assumptions; coarse temporal resolution

---

### ConnectivityKernel

**Input:** k-nearest neighbor graph
**Biological signal:** Transcriptomic similarity
**When to use:** As a smoothing kernel combined with directional kernels

```python
from cellrank.kernels import ConnectivityKernel

ck = ConnectivityKernel(adata)
ck.compute_transition_matrix()
```

**Note:** ConnectivityKernel alone has no directionality. Use it to smooth noisy directional kernels:
```python
combined = 0.8 * directional_kernel + 0.2 * ck
```

---

## Kernel Combination

### Syntax

```python
# Two kernels
combined = 0.8 * vk + 0.2 * ck

# Three kernels
combined = 0.5 * vk + 0.3 * ctk + 0.2 * ck

# Equal weights
combined = 0.5 * vk + 0.5 * ctk
```

### Weight Selection Guidelines

| Scenario | Suggested Weights |
|----------|-------------------|
| High-quality velocity | `0.8 * vk + 0.2 * ck` |
| Moderate velocity quality | `0.6 * vk + 0.2 * ctk + 0.2 * ck` |
| Poor velocity, good CytoTRACE | `0.3 * vk + 0.5 * ctk + 0.2 * ck` |
| No velocity | `0.8 * ctk + 0.2 * ck` |
| Time-course + velocity | `0.5 * rtk + 0.3 * vk + 0.2 * ck` |
| Pseudotime + CytoTRACE | `0.5 * ptk + 0.3 * ctk + 0.2 * ck` |

### Best Practice
1. Run each directional kernel individually first
2. Compare terminal states and fate probabilities
3. Combine with weights that reflect data quality
4. Report both individual and combined results

---

## GPCCA Estimator

The **Generalized Perron Cluster Cluster Analysis** estimator is used for all downstream analysis, regardless of which kernel(s) you use.

### Full Pipeline

```python
from cellrank.estimators import GPCCA

# 1. Create estimator from any kernel
g = GPCCA(kernel)

# 2. Schur decomposition — eigenvalue analysis
g.compute_schur(n_components=20)
g.plot_spectrum(real_only=True)  # Look for eigenvalue gap

# 3. Macrostates — metastable cell populations
g.compute_macrostates(n_states=6, cluster_key='clusters')
g.plot_macrostates(which='all', basis='umap')

# 4. Terminal states — absorbing endpoints
g.set_terminal_states()  # Automatic selection
# OR: g.set_terminal_states(states=['State_A', 'State_B'])  # Manual
g.plot_macrostates(which='terminal', basis='umap')

# 5. Fate probabilities — absorption probabilities
g.compute_fate_probabilities()
g.plot_fate_probabilities(basis='umap')

# 6. Driver genes — lineage-specific regulators
drivers = g.compute_lineage_drivers(lineages='State_A', return_drivers=True)
print(drivers.head(20))

# 7. Gene trends — expression along fate
model = cr.models.GAM(adata)
cr.pl.gene_trends(adata, model=model, genes=['Gene1', 'Gene2'])
```

### Choosing n_states

- Plot the eigenvalue spectrum with `g.plot_spectrum(real_only=True)`
- Look for a **gap** in the real eigenvalues
- The number of eigenvalues close to 1 (before the gap) suggests the number of metastable macrostates
- Try a few values and check biological interpretability

---

## Decision Flowchart

```
START: What data do I have?
│
├── Spliced + unspliced counts?
│   ├── YES, good velocity quality
│   │   └── Primary: VelocityKernel + ConnectivityKernel
│   │       Consider adding: CytoTRACEKernel for validation
│   │
│   └── YES, but poor velocity quality
│       └── Primary: CytoTRACEKernel
│           Consider: PseudotimeKernel as alternative
│           Optional: VelocityKernel with low weight
│
├── Only standard scRNA-seq (no velocity layers)?
│   └── Primary: CytoTRACEKernel
│       Alternative: PseudotimeKernel (from DPT or other method)
│
├── Time-course experiment?
│   └── Primary: RealTimeKernel
│       Combine with: VelocityKernel or CytoTRACEKernel if available
│
└── Metabolic labeling data?
    └── Process with dynamo/scVelo → VelocityKernel
        Combine with: CytoTRACEKernel for validation
```

---

## Common Issues & Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| No terminal states identified | Too few macrostates or wrong `n_states` | Adjust `n_states`; check eigenvalue gap |
| All cells assigned to one fate | Data may not have branching | Check if trajectory is linear; try fewer macrostates |
| GPCCA fails with error | Transition matrix issues | Check kernel computation; try `softmax_scale` parameter |
| Kernels disagree strongly | Different data views see different biology | Report both; investigate biologically |
| Very slow computation | Large cell numbers | Subsample or use `n_components` < 20 |
| Driver genes don't make sense | Wrong terminal states or noisy data | Validate terminal states with markers; try combined kernel |

---

## Citing CellRank

When using CellRank in publications, cite:

1. **CellRank 1:** Lange et al. (2022) "CellRank for directed single-cell fate mapping." *Nat. Methods* 19:159-170. [doi:10.1038/s41592-021-01346-6](https://doi.org/10.1038/s41592-021-01346-6)

2. **CellRank 2:** Weiler et al. (2024) "CellRank 2: unified fate mapping in multiview single-cell data." *Nat. Methods* 21:1196-1205. [doi:10.1038/s41592-024-02303-9](https://doi.org/10.1038/s41592-024-02303-9)

3. **Protocol:** Weiler & Theis (2026) "CellRank: consistent and data view agnostic fate mapping." *Nat. Protocols*. [doi:10.1038/s41596-025-01314-w](https://doi.org/10.1038/s41596-025-01314-w)

4. **pyGPCCA:** Reuter et al. (2019) *J. Chem. Phys.* 150:174103

---

## External Resources

| Resource | URL |
|----------|-----|
| CellRank documentation | [cellrank.readthedocs.io](https://cellrank.readthedocs.io/) |
| CellRank GitHub | [github.com/theislab/cellrank](https://github.com/theislab/cellrank) |
| CellRank Protocol repo | [github.com/theislab/cellrank_protocol](https://github.com/theislab/cellrank_protocol) |
| Protocol Jupyter Book | [theislab.github.io/cellrank_protocol](https://theislab.github.io/cellrank_protocol/index.html) |
| Protocol data (Figshare) | [doi.org/10.6084/m9.figshare.c.7752290.v1](https://doi.org/10.6084/m9.figshare.c.7752290.v1) |
| scverse forum | [discourse.scverse.org/c/ecosystem/cellrank](https://discourse.scverse.org/c/ecosystem/cellrank/) |
| Best practices paper | Heumos et al. (2023) *Nat. Rev. Genet.* 24:550-572 |
