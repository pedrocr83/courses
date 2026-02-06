# Course: Trajectory & Pseudotime Analysis

## Understanding Cellular Dynamics and Fate Mapping

---

## Course Overview

| | |
|---|---|
| **Duration** | 4 weeks (20-25 hours total) |
| **Format** | Self-paced with hands-on labs |
| **Level** | Intermediate to Advanced |
| **Prerequisites** | Clustering & annotation, basic statistics |
| **Updated** | February 2026 — includes CellRank 2 multiview framework |

## Start Here (Do This First)

- **Step-by-step checklist**: `START_HERE.md`
- **Glossary**: `glossary.md`
- **Environment setup**: `setup.md`

### Learning Objectives

By the end of this course, you will be able to:
1. Understand when trajectory analysis is appropriate
2. Infer pseudotime ordering using multiple methods
3. Construct and interpret lineage trees
4. Identify genes that change along trajectories
5. Apply RNA velocity for directional inference
6. **Use CellRank 2 for multiview fate mapping** — combining pseudotime, RNA velocity, developmental potential, and experimental time points
7. **Compare CellRank kernels** (VelocityKernel, CytoTRACEKernel, PseudotimeKernel, RealTimeKernel) and choose the appropriate one for your data
8. **Perform kernel combination** to leverage complementary data views
9. **Infer terminal/initial states, fate probabilities, and driver genes** in a scalable, reproducible manner

---

## Course Structure

### Week 1: Foundations of Trajectory Analysis

#### Module 1: When Cells Form Trajectories (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 1.1 | Continuous vs discrete cell states | Lecture |
| 1.2 | Biological contexts: development, differentiation, activation | Lecture |
| 1.3 | Snapshot data and temporal inference | Lecture |
| **Lab 1** | Explore developmental scRNA-seq data | Hands-on |

**Key Concepts:** Manifold, continuous processes, temporal ordering

---

#### Module 2: Pseudotime Concepts (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 2.1 | What is pseudotime? | Lecture |
| 2.2 | Assumptions and limitations | Lecture |
| 2.3 | Root cell selection | Lecture |
| 2.4 | Linear vs branching trajectories | Lecture |
| **Lab 2** | Visualize pseudotime on known developmental data | Hands-on |

**Key Concepts:** Ordering cells, directionality, branching points

---

#### Module 3: Diffusion-Based Methods (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 3.1 | Diffusion maps and diffusion pseudotime | Lecture |
| 3.2 | PAGA: partition-based graph abstraction | Lecture |
| 3.3 | Combining PAGA with embeddings | Demo |
| **Lab 3** | Run diffusion pseudotime with Scanpy | Hands-on |

**Key Concepts:** Transition probabilities, diffusion components, PAGA graphs

---

### Week 2: Trajectory Inference Methods

#### Module 4: Monocle 3 (2.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 4.1 | UMAP + graph learning | Lecture |
| 4.2 | Principal graph and pseudotime | Lecture |
| 4.3 | Choosing root nodes | Demo |
| 4.4 | Multiple lineages and branches | Lecture |
| **Lab 4** | Complete Monocle 3 workflow | Hands-on |

**Key Concepts:** Learn_graph, order_cells, branches

---

#### Module 5: Slingshot (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 5.1 | Cluster-based trajectory inference | Lecture |
| 5.2 | Lineage identification | Lecture |
| 5.3 | Smooth curves through clusters | Lecture |
| **Lab 5** | Run Slingshot on developmental data | Hands-on |

**Key Concepts:** Minimum spanning tree, principal curves, lineages

---

#### Module 6: RNA Velocity (2.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 6.1 | Spliced vs unspliced RNA | Lecture |
| 6.2 | Velocity estimation | Lecture |
| 6.3 | scVelo: stochastic and dynamical models | Lecture |
| 6.4 | Interpreting velocity arrows | Demo |
| 6.5 | Limitations & current challenges in RNA velocity | Lecture |
| **Lab 6** | Compute RNA velocity with scVelo | Hands-on |

**Key Concepts:** Splicing kinetics, future state prediction, velocity embedding

**Key References:**
- La Manno et al. (2018) — RNA velocity of single cells, *Nature*
- Bergen et al. (2020) — Generalizing RNA velocity, *Nat. Biotechnol.*
- Bergen et al. (2021) — RNA velocity: current challenges, *Mol. Syst. Biol.*

---

### Week 3: CellRank 2 — Multiview Fate Mapping

> **Based on:** Weiler & Theis (2026) "CellRank: consistent and data view agnostic fate mapping for single-cell genomics", *Nature Protocols*. [doi:10.1038/s41596-025-01314-w](https://doi.org/10.1038/s41596-025-01314-w)
>
> **Code:** [github.com/theislab/cellrank](https://github.com/theislab/cellrank) | [github.com/theislab/cellrank_protocol](https://github.com/theislab/cellrank_protocol)

#### Module 7: CellRank 2 Framework Overview (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 7.1 | From CellRank 1 to CellRank 2: evolution of the framework | Lecture |
| 7.2 | The kernel abstraction: transition matrices from diverse data views | Lecture |
| 7.3 | Markov chains and GPCCA for macrostate identification | Lecture |
| 7.4 | CellRank architecture: Kernels → Estimators → Analysis | Demo |
| **Lab 7** | CellRank 2 overview — load data & explore the API | Hands-on |

**Key Concepts:**
- **Kernels**: Objects that compute a cell-cell transition matrix from a specific data view
- **Estimators**: GPCCA-based analysis of the transition matrix (macrostates, fate probabilities)
- **Modularity**: Any kernel can be combined with any estimator
- **pyGPCCA**: Backend for Generalized Perron Cluster Cluster Analysis

---

#### Module 8: VelocityKernel — RNA Velocity-Based Fate Mapping (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 8.1 | Building a transition matrix from RNA velocity | Lecture |
| 8.2 | Combining VelocityKernel with ConnectivityKernel | Lecture |
| 8.3 | Terminal state identification via GPCCA | Demo |
| 8.4 | Fate probabilities, absorption times, and driver genes | Lecture |
| **Lab 8** | CellRank VelocityKernel on pancreas data | Hands-on |

**Key Concepts:** Velocity-based transitions, kernel weighting, Schur decomposition, macrostates

**Lab Reference:** Corresponds to `notebooks/velocity/` in the [CellRank Protocol](https://theislab.github.io/cellrank_protocol/)

---

#### Module 9: CytoTRACEKernel — Developmental Potential (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 9.1 | CytoTRACE: estimating developmental potential from gene counts | Lecture |
| 9.2 | CytoTRACEKernel: building transitions from stemness scores | Lecture |
| 9.3 | When to use CytoTRACEKernel vs other kernels | Discussion |
| **Lab 9A** | CytoTRACEKernel on bone marrow data | Hands-on |

**Key Concepts:** Gene count complexity as stemness proxy, directionality without velocity

**Lab Reference:** Corresponds to `notebooks/cytotrace/` in the CellRank Protocol

---

#### Module 10: PseudotimeKernel — Pseudotime-Driven Fate Mapping (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 10.1 | Using any pseudotime as input to CellRank | Lecture |
| 10.2 | Automatic vs manual terminal state definition | Lecture |
| 10.3 | Comparing CytoTRACEKernel and PseudotimeKernel results | Demo |
| **Lab 9B** | PseudotimeKernel: DPT → CellRank fate probabilities | Hands-on |

**Key Concepts:** Pseudotime as directional prior, kernel comparison, reproducibility

**Lab Reference:** Corresponds to `notebooks/pseudotime/` in the CellRank Protocol

---

#### Module 11: RealTimeKernel — Experimental Time Points (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 11.1 | When you have real experimental time: time-course scRNA-seq | Lecture |
| 11.2 | Optimal transport for time-resolved data | Lecture |
| 11.3 | RealTimeKernel: connecting cells across measured time points | Lecture |
| 11.4 | Metabolic labeling as a special case | Lecture |
| **Lab 9C** | RealTimeKernel on time-course data | Hands-on |

**Key Concepts:** Optimal transport, Waddington-OT, time-resolved fate mapping, metabolic labeling (scSLAM-seq, sci-fate)

**Lab Reference:** Corresponds to `notebooks/realtime/` and `notebooks/metabolic_labeling/` in the CellRank Protocol

---

#### Module 12: Kernel Combination & Multiview Analysis (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 12.1 | Why combine kernels? Complementary views of the same biology | Lecture |
| 12.2 | Weighted kernel addition: `0.8 * vk + 0.2 * ck` | Demo |
| 12.3 | Evaluating kernel agreement and disagreement | Lecture |
| 12.4 | Full multiview workflow: velocity + pseudotime + CytoTRACE | Lecture |
| **Lab 9D** | Kernel combination: compare single vs combined kernel results | Hands-on |

**Key Concepts:** Kernel arithmetic, multiview integration, view-agnostic fate mapping, consistency checks

**Lab Reference:** Corresponds to the combined analysis in the CellRank Protocol

---

### Week 4: Downstream Analysis & Advanced Topics

#### Module 13: Trajectory Differential Expression (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 13.1 | Genes changing along pseudotime | Lecture |
| 13.2 | GAMs and spline regression | Lecture |
| 13.3 | tradeSeq for trajectory DE | Lecture |
| 13.4 | CellRank driver genes: lineage-specific regulators | Lecture |
| **Lab 10** | Find pseudotime-associated genes + CellRank driver genes | Hands-on |

**Key Concepts:** Smooth regression, branch-specific DE, gene modules, driver gene analysis

---

#### Module 14: Gene Expression Trends & Fate Visualization (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 14.1 | Expression heatmaps along pseudotime | Lecture |
| 14.2 | CellRank gene trends: expression along fate probabilities | Lecture |
| 14.3 | Clustering gene trends for module discovery | Demo |
| **Lab 11** | Visualize gene dynamics: classical + CellRank trend plots | Hands-on |

**Key Concepts:** Ordered heatmaps, gene trends along lineages, gene module clustering

---

#### Module 15: Method Comparison & Best Practices (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 15.1 | Benchmark studies overview (dynverse, CellRank protocol comparisons) | Reading |
| 15.2 | Choosing the right method and kernel for your data | Discussion |
| 15.3 | Validation strategies: biological, computational, experimental | Lecture |
| 15.4 | Reporting standards for trajectory & fate mapping results | Lecture |
| **Lab 12** | Compare methods and kernels on same dataset | Hands-on |

**Key Concepts:** dynverse benchmark, CellRank kernel comparison, ground truth, biological validation

---

## Assessments

### Quizzes (after key modules)
- 5-10 questions covering key concepts
- Immediate feedback

| Quiz | After Module | Topic |
|------|-------------|-------|
| Quiz 1 | Module 2 | Pseudotime concepts & assumptions |
| Quiz 2 | Module 6 | RNA velocity theory & interpretation |
| Quiz 3 | Module 7-8 | CellRank 2 framework & VelocityKernel |
| Quiz 4 | Module 12 | Multiview kernels & kernel combination |

### Practical Assignments

| Assignment | Module | Description |
|------------|--------|-------------|
| **A1** | 3-4 | Pseudotime inference comparison |
| **A2** | 6 | RNA velocity analysis |
| **A3** | 13-14 | Trajectory DE and visualization |
| **A4** | 9-12 | **CellRank 2 multiview fate mapping** |

### Final Project
**Complete trajectory and fate mapping analysis of differentiation data**
- Infer pseudotime with 2+ methods
- Compute RNA velocity
- Apply CellRank 2 with at least 2 different kernels
- Compare single-kernel vs combined-kernel fate probabilities
- Identify terminal states, fate probabilities, and driver genes
- Biological interpretation and validation strategy

---

## Software Setup

### Python (Primary: CellRank Protocol Environment)

```bash
# Recommended: Use the CellRank Protocol environment
conda create -n cellrank-course python=3.11 --yes
conda activate cellrank-course
conda install -c conda-forge cellrank

# Clone the protocol repo for reference notebooks
git clone https://github.com/theislab/cellrank_protocol.git
cd cellrank_protocol
pip install -e ".[jupyter]"
cd ..

# Additional tools
pip install scvelo loompy
pip install jupyterlab

# Register kernel
python -m ipykernel install --user --name cellrank-course --display-name "cellrank-course"
```

### Verify Installation

```python
import scanpy as sc
import scvelo as scv
import cellrank as cr

print(f"scanpy:   {sc.__version__}")
print(f"scvelo:   {scv.__version__}")
print(f"cellrank: {cr.__version__}")

# CellRank 2 kernels should be available
from cellrank.kernels import VelocityKernel, CytoTRACEKernel, PseudotimeKernel, RealTimeKernel
print("All CellRank 2 kernels available!")
```

### R (Optional: Monocle3, Slingshot)

```r
# Monocle3
devtools::install_github("cole-trapnell-lab/monocle3")

# Slingshot
BiocManager::install("slingshot")

# tradeSeq for DE
BiocManager::install("tradeSeq")
```

---

## Demo Datasets

| Dataset | Process | Cells | Use | CellRank Kernel |
|---------|---------|-------|-----|-----------------|
| Pancreas endocrine | Differentiation | 3k | Primary labs | VelocityKernel |
| Human bone marrow | Multi-lineage hematopoiesis | 5k | CytoTRACE + Pseudotime labs | CytoTRACEKernel, PseudotimeKernel |
| Time-course reprogramming | Real-time trajectory | Variable | RealTimeKernel lab | RealTimeKernel |
| Cell cycle | Cyclic | 1k | RNA velocity exploration | VelocityKernel |

### Data Access

```python
# Protocol datasets (Figshare)
# https://doi.org/10.6084/m9.figshare.c.7752290.v1

# Built-in scVelo datasets
import scvelo as scv
adata_pancreas = scv.datasets.pancreas()
adata_bonemarrow = scv.datasets.bonemarrow()
adata_dentategyrus = scv.datasets.dentategyrus()
```

---

## Methods Summary

| Method | Approach | Branches | Direction | CellRank Kernel |
|--------|----------|----------|-----------|-----------------|
| Diffusion PT | Random walks | No | Needs root | Input to PseudotimeKernel |
| PAGA | Graph abstraction | Yes | Needs root | — |
| Monocle 3 | Principal graph | Yes | Needs root | — |
| Slingshot | MST + curves | Yes | Needs root | — |
| RNA Velocity | Splicing dynamics | Implicit | Intrinsic | VelocityKernel |
| CytoTRACE | Gene count complexity | No | Intrinsic | CytoTRACEKernel |
| Experimental time | Sample time points | No | Known | RealTimeKernel |
| **CellRank 2** | **Markov chains (any kernel)** | **Yes** | **From kernel** | **Any / Combined** |

---

## CellRank 2 Kernel Decision Guide

```
Which data views do you have?
│
├── RNA velocity (spliced/unspliced counts)?
│   └── YES → VelocityKernel
│       └── Consider combining with ConnectivityKernel
│
├── Developmental potential / stemness information?
│   └── YES → CytoTRACEKernel
│       └── No velocity needed; uses gene count complexity
│
├── Pseudotime from any method (DPT, Monocle, Slingshot)?
│   └── YES → PseudotimeKernel
│       └── Wraps any pseudotime into CellRank framework
│
├── Experimental time points (time-course experiment)?
│   └── YES → RealTimeKernel
│       └── Optimal transport between time points
│
└── Multiple views available?
    └── YES → Combine kernels!
        └── Example: 0.8 * VelocityKernel + 0.2 * CytoTRACEKernel
        └── Evaluate each kernel individually first, then combine
```

---

## Common Mistakes to Avoid

| Mistake | Consequence | Prevention |
|---------|-------------|------------|
| Trajectory on discrete states | Meaningless pseudotime | Verify biological continuity |
| Wrong root selection | Reversed trajectory | Use biological knowledge |
| Ignoring batch effects | Technical pseudotime | Integrate first |
| Over-interpreting velocity | False directionality | Validate with known markers |
| Single method | Biased conclusions | Compare multiple approaches |
| No gene validation | Unvalidated ordering | Check known temporal genes |
| Using only VelocityKernel when velocity quality is poor | Unreliable fate mapping | Check scVelo diagnostics; try CytoTRACEKernel or PseudotimeKernel |
| Not comparing kernels | Missing complementary signals | Run ≥2 kernels and compare fate probabilities |
| Ignoring CellRank's GPCCA macrostate analysis | Ad-hoc terminal state selection | Use `compute_macrostates()` for data-driven states |

---

## What Bad Trajectory Analysis Looks Like (and How to Catch It)

### Symptoms (Red Flags)

- **No continuity**: pseudotime jumps between unrelated clusters; ordering contradicts known biology.
- **Root sensitivity**: small changes in root selection reverse the story.
- **Batch-driven trajectory**: pseudotime correlates with batch/donor rather than known stages.
- **Velocity chaos**: velocity arrows point everywhere; no consistent flow; velocity genes are unstable.
- **Kernel disagreement**: CytoTRACEKernel and VelocityKernel give opposite fate probabilities with no biological explanation.

### Common Failure Modes → What It Looks Like → Fix

- **Applying trajectory to discrete cell types**
  - **Looks like**: "trajectory" is just cluster ordering with no transitional cells.
  - **Fix**: restrict to a lineage/compartment; re-check whether biology is continuous.

- **Wrong root / terminal definitions**
  - **Looks like**: late markers appear early; known progenitors show high pseudotime.
  - **Fix**: choose root using prior knowledge + marker genes; test multiple plausible roots and report sensitivity.

- **Batch effects not addressed**
  - **Looks like**: pseudotime aligns with batches; stages are batch-specific.
  - **Fix**: integrate/correct batch before trajectory (or analyze per-batch and compare).

- **Over-interpreting RNA velocity**
  - **Looks like**: arrows interpreted as "fate" without checking spliced/unspliced quality or model fit.
  - **Fix**: use scVelo diagnostics; validate direction using known stage markers; treat velocity as hypothesis, not proof.

- **CellRank-specific: Wrong number of macrostates**
  - **Looks like**: GPCCA merges distinct terminal populations, or splits one population artificially.
  - **Fix**: Inspect Schur decomposition eigenvalue gap; try different `n_states`; validate macrostates against known markers.

### Minimum "Sanity Checks" Before You Report

- **Biology**: known early/late markers change in the expected direction along pseudotime.
- **Robustness**: major conclusions hold across (a) method choice, (b) kernel choice, and (c) reasonable parameter changes.
- **Confounding**: pseudotime is not primarily explained by batch/donor.
- **Kernel consistency**: if using CellRank 2, compare at least 2 kernels — agreement strengthens conclusions.
- **Validation plan**: you can propose an orthogonal validation (perturbation, lineage tracing, protein).

---

## Learning Path Progression

```
BEGINNER (Week 1)
│
├── Understand trajectory concepts
├── Run diffusion pseudotime
└── Visualize pseudotime on UMAP
        │
        ▼
INTERMEDIATE (Week 2)
│
├── Monocle 3 / Slingshot
├── Handle branching trajectories
├── RNA velocity basics
└── Interpret velocity plots
        │
        ▼
ADVANCED (Week 3: CellRank 2)
│
├── CellRank 2 framework & kernels
├── VelocityKernel fate mapping
├── CytoTRACEKernel developmental potential
├── PseudotimeKernel & RealTimeKernel
└── Kernel combination & multiview analysis
        │
        ▼
EXPERT (Week 4+)
│
├── Trajectory differential expression
├── CellRank driver genes & gene trends
├── Multi-kernel comparison & validation
└── Publication-ready fate mapping
```

---

## Connection to Other Courses

| Prerequisite from | Why |
|-------------------|-----|
| Clustering & Annotation | Need cell type context for terminal state interpretation |
| Data Integration | For multi-batch trajectory; batch correction before fate mapping |

| Enables | Why |
|---------|-----|
| Differential Expression | DE along pseudotime; CellRank driver genes |
| Cell-Cell Communication | Signaling dynamics during differentiation |

---

## Biological Applications

- **Development:** Organogenesis, cell fate decisions, lineage commitment
- **Stem cells:** Differentiation pathways, reprogramming dynamics
- **Cancer:** Tumor evolution, EMT trajectories, drug resistance
- **Immunology:** T cell activation, hematopoietic differentiation
- **Regeneration:** Tissue repair dynamics, wound healing
- **Disease progression:** Neurodegeneration, fibrosis temporal dynamics

---

## Key References

### Foundational
- La Manno et al. (2018) "RNA velocity of single cells" *Nature* 560:494-498
- Bergen et al. (2020) "Generalizing RNA velocity to transient cell states" *Nat. Biotechnol.* 38:1408-1414
- Haghverdi et al. (2016) "Diffusion pseudotime robustly reconstructs lineage branching" *Nat. Methods* 13:845-848

### CellRank
- Lange et al. (2022) "CellRank for directed single-cell fate mapping" *Nat. Methods* 19:159-170
- Weiler et al. (2024) "CellRank 2: unified fate mapping in multiview single-cell data" *Nat. Methods* 21:1196-1205
- **Weiler & Theis (2026) "CellRank: consistent and data view agnostic fate mapping" *Nat. Protocols*** [doi:10.1038/s41596-025-01314-w](https://doi.org/10.1038/s41596-025-01314-w)

### Best Practices
- Heumos et al. (2023) "Best practices for single-cell analysis across modalities" *Nat. Rev. Genet.* 24:550-572
- Weiler et al. (2023) "A guide to trajectory inference and RNA velocity" *Methods Mol. Biol.* 2584:269-292
- Bergen et al. (2021) "RNA velocity — current challenges and future perspectives" *Mol. Syst. Biol.* 17:e10282

### Protocol & Data
- CellRank Protocol: [github.com/theislab/cellrank_protocol](https://github.com/theislab/cellrank_protocol)
- Protocol Jupyter Book: [theislab.github.io/cellrank_protocol](https://theislab.github.io/cellrank_protocol/index.html)
- Protocol data: [Figshare collection](https://doi.org/10.6084/m9.figshare.c.7752290.v1)

---

*This course reveals the temporal dimension hidden in snapshot single-cell data — now with CellRank 2's unified, multiview approach to fate mapping.*
