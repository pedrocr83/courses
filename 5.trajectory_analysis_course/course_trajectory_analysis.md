# Course: Trajectory & Pseudotime Analysis

## Understanding Cellular Dynamics and Differentiation

---

## Course Overview

| | |
|---|---|
| **Duration** | 3 weeks (15-20 hours total) |
| **Format** | Self-paced with hands-on labs |
| **Level** | Intermediate to Advanced |
| **Prerequisites** | Clustering & annotation, basic statistics |

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
| **Lab 6** | Compute RNA velocity with scVelo | Hands-on |

**Key Concepts:** Splicing kinetics, future state prediction, velocity embedding

---

### Week 3: Downstream Analysis & Advanced Topics

#### Module 7: Trajectory Differential Expression (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 7.1 | Genes changing along pseudotime | Lecture |
| 7.2 | GAMs and spline regression | Lecture |
| 7.3 | tradeSeq for trajectory DE | Lecture |
| **Lab 7** | Find pseudotime-associated genes | Hands-on |

**Key Concepts:** Smooth regression, branch-specific DE, gene modules

---

#### Module 8: Gene Dynamics Visualization (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 8.1 | Expression heatmaps along pseudotime | Lecture |
| 8.2 | Gene cascades and modules | Lecture |
| 8.3 | Branch-specific expression | Demo |
| **Lab 8** | Visualize gene dynamics | Hands-on |

**Key Concepts:** Ordered heatmaps, gene trends, branch comparison

---

#### Module 9: CellRank: Fate Probabilities (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 9.1 | Markov chains for fate prediction | Lecture |
| 9.2 | Combining velocity with transcriptomics | Lecture |
| 9.3 | Terminal states and fate probabilities | Lecture |
| **Lab 9** | Predict cell fates with CellRank | Hands-on |

**Key Concepts:** Transition matrix, absorption probabilities, fate bias

---

#### Module 10: Method Comparison & Best Practices (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 10.1 | Benchmark studies overview | Reading |
| 10.2 | Choosing the right method | Discussion |
| 10.3 | Validation strategies | Lecture |
| **Lab 10** | Compare methods on same dataset | Hands-on |

**Key Concepts:** dynverse benchmark, ground truth, biological validation

---

## Assessments

### Quizzes (after each module)
- 5-10 questions covering key concepts
- Immediate feedback

### Practical Assignments

| Assignment | Module | Description |
|------------|--------|-------------|
| **A1** | 3-4 | Pseudotime inference comparison |
| **A2** | 6 | RNA velocity analysis |
| **A3** | 7-8 | Trajectory DE and visualization |

### Final Project
**Complete trajectory analysis of differentiation data**
- Infer pseudotime with 2+ methods
- Compute RNA velocity
- Identify trajectory-associated genes
- Biological interpretation

---

## Software Setup

### Python (Scanpy, scVelo)

```python
pip install scanpy scvelo cellrank
pip install loompy  # For velocity data
```

### R (Monocle3, Slingshot)

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

| Dataset | Process | Cells | Use |
|---------|---------|-------|-----|
| Pancreas endocrine | Differentiation | 3k | Primary labs |
| Hematopoiesis | Multi-lineage | 5k | Branching |
| Cell cycle | Cyclic | 1k | RNA velocity |

---

## Methods Summary

| Method | Approach | Branches | Direction |
|--------|----------|----------|-----------|
| Diffusion PT | Random walks | No | Needs root |
| PAGA | Graph abstraction | Yes | Needs root |
| Monocle 3 | Principal graph | Yes | Needs root |
| Slingshot | MST + curves | Yes | Needs root |
| RNA Velocity | Splicing dynamics | Implicit | Intrinsic |
| CellRank | Markov chains | Yes | From velocity |

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

---

## What Bad Trajectory Analysis Looks Like (and How to Catch It)

### Symptoms (Red Flags)

- **No continuity**: pseudotime jumps between unrelated clusters; ordering contradicts known biology.
- **Root sensitivity**: small changes in root selection reverse the story.
- **Batch-driven trajectory**: pseudotime correlates with batch/donor rather than known stages.
- **Velocity chaos**: velocity arrows point everywhere; no consistent flow; velocity genes are unstable.

### Common Failure Modes → What It Looks Like → Fix

- **Applying trajectory to discrete cell types**  
  - **Looks like**: “trajectory” is just cluster ordering with no transitional cells.  
  - **Fix**: restrict to a lineage/compartment; re-check whether biology is continuous.

- **Wrong root / terminal definitions**  
  - **Looks like**: late markers appear early; known progenitors show high pseudotime.  
  - **Fix**: choose root using prior knowledge + marker genes; test multiple plausible roots and report sensitivity.

- **Batch effects not addressed**  
  - **Looks like**: pseudotime aligns with batches; stages are batch-specific.  
  - **Fix**: integrate/correct batch before trajectory (or analyze per-batch and compare).

- **Over-interpreting RNA velocity**  
  - **Looks like**: arrows interpreted as “fate” without checking spliced/unspliced quality or model fit.  
  - **Fix**: use scVelo diagnostics; validate direction using known stage markers; treat velocity as hypothesis, not proof.

### Minimum “Sanity Checks” Before You Report

- **Biology**: known early/late markers change in the expected direction along pseudotime.
- **Robustness**: major conclusions hold across (a) method choice and (b) reasonable parameter changes.
- **Confounding**: pseudotime is not primarily explained by batch/donor.
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
ADVANCED (Week 3+)
│
├── Trajectory differential expression
├── CellRank fate probabilities
├── Multi-modal trajectory analysis
└── Custom method development
```

---

## Connection to Other Courses

| Prerequisite from | Why |
|-------------------|-----|
| Clustering & Annotation | Need cell type context |
| Data Integration | For multi-batch trajectory |

| Enables | Why |
|---------|-----|
| Differential Expression | DE along pseudotime |
| Cell-Cell Communication | Signaling during differentiation |

---

## Biological Applications

- **Development:** Organogenesis, cell fate decisions
- **Stem cells:** Differentiation pathways
- **Cancer:** Tumor evolution, EMT
- **Immunology:** T cell activation, differentiation
- **Regeneration:** Tissue repair dynamics

---

*This course reveals the temporal dimension hidden in snapshot single-cell data.*

