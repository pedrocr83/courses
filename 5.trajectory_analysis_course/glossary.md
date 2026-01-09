# Glossary: Trajectory & Pseudotime Analysis

---

## Trajectory Concepts

- **Cell state**: A transcriptional configuration (often continuous rather than discrete).
- **Trajectory**: Continuous progression of states through a biological process.
- **Pseudotime**: A **relative** ordering of cells along a trajectory (not real time).
- **Root cell**: Starting point for pseudotime ordering (must be chosen/estimated).
- **Terminal state**: End state(s) of a process; may branch into multiple fates.
- **Branching**: Multiple lineage outcomes from a shared progenitor-like state.

---

## Methods

- **Diffusion map**: Embedding using random-walk transition probabilities; captures manifold structure.
- **DPT (Diffusion Pseudotime)**: Pseudotime computed using diffusion distances from a root.
- **PAGA**: Graph abstraction summarizing connectivity between clusters; supports branching structure.
- **Monocle 3**: Learns a principal graph in low-dimensional space and orders cells along it.
- **Slingshot**: Fits smooth curves (principal curves) through cluster-based lineage trees.

---

## RNA Velocity

- **Spliced / unspliced counts**: Mature mRNA vs pre-mRNA; used to infer transcriptional dynamics.
- **RNA velocity**: Predicts the **direction** of expression change (a short-term “future state”).
- **Velocity graph**: Graph encoding transition likelihoods informed by velocity vectors.
- **Stochastic model (scVelo)**: Simpler velocity model; faster, fewer kinetics assumptions.
- **Dynamical model (scVelo)**: Fits gene-specific kinetics; more expressive but more demanding.

---

## Fate Inference

- **Markov chain**: Probabilistic model of transitions between states.
- **Transition matrix**: Cells × cells probabilities of moving from one state to another.
- **CellRank**: Combines velocity and transcriptomic similarity to compute fate probabilities.
- **Macrostates**: Metastable coarse-grained states inferred by CellRank.
- **Fate probability**: Probability a cell ends in each terminal state.
- **Driver genes**: Genes correlated with fate probabilities (candidate regulators).

---

## Where You’ll Use These

- **DPT/diffusion map/root selection**: `labs/lab03_diffusion_pt.ipynb`
- **Velocity, spliced/unspliced**: `labs/lab06_velocity.ipynb`
- **Trajectory-associated genes**: `labs/lab07_trajectory_de.ipynb`
- **CellRank/macrostates/fates**: `labs/lab09_cellrank.ipynb`
- **Method comparison**: `labs/lab10_comparison.ipynb`


