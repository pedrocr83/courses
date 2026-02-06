# Quiz: CellRank 2 — Multiview Fate Mapping

**Covers Modules 7-12** — CellRank 2 framework, kernels, and kernel combination

---

## Part A: CellRank 2 Framework (Modules 7-8)

### Q1. What is a "kernel" in CellRank 2?
- a) A machine learning algorithm for classification
- b) An object that computes a cell-cell transition matrix from a specific data view
- c) A gene expression clustering method
- d) A Python package for deep learning

<details><summary>Answer</summary>b) A kernel in CellRank 2 is an object that computes a cell-cell transition matrix from a specific data view (e.g., RNA velocity, pseudotime, developmental potential).</details>

---

### Q2. What does GPCCA stand for, and what does it do?
- a) General Purpose Cell Cluster Analysis — clusters cells
- b) Generalized Perron Cluster Cluster Analysis — identifies metastable macrostates from the transition matrix
- c) Gene-based Principal Component Correlation Analysis — reduces dimensionality
- d) Graph-Partitioned Cell Classification Algorithm — classifies cell types

<details><summary>Answer</summary>b) Generalized Perron Cluster Cluster Analysis — it identifies metastable macrostates (groups of cells the Markov chain tends to stay in) from the transition matrix computed by a kernel.</details>

---

### Q3. In CellRank 2, what are "fate probabilities"?
- a) The probability that a gene is expressed in a cell
- b) The absorption probability of each cell reaching each terminal state
- c) The likelihood that a cell type exists in the sample
- d) The p-value from differential expression testing

<details><summary>Answer</summary>b) Fate probabilities are the absorption probabilities — for each cell, the probability of eventually reaching each terminal (absorbing) state in the Markov chain. They sum to 1 across all terminal states.</details>

---

### Q4. What is the recommended approach when combining VelocityKernel with ConnectivityKernel?
- a) Always use equal weights (0.5 + 0.5)
- b) Give higher weight to the directional kernel (e.g., 0.8 velocity + 0.2 connectivity)
- c) Only use ConnectivityKernel; VelocityKernel is not needed
- d) Multiply them together: vk * ck

<details><summary>Answer</summary>b) Give higher weight to the directional kernel (VelocityKernel) to preserve the biological signal, while using ConnectivityKernel at lower weight to smooth noise. Typical: 0.8 * vk + 0.2 * ck.</details>

---

### Q5. How do you determine the appropriate number of macrostates in GPCCA?
- a) Always use n_states=2
- b) Use the number of clusters from Leiden clustering
- c) Examine the eigenvalue spectrum gap from the Schur decomposition
- d) Count the number of genes in the dataset

<details><summary>Answer</summary>c) Plot the eigenvalue spectrum with `g.plot_spectrum(real_only=True)` and look for a gap in the real eigenvalues close to 1. The number of eigenvalues before the gap suggests the number of metastable macrostates.</details>

---

## Part B: Individual Kernels (Modules 9-11)

### Q6. The CytoTRACEKernel infers developmental direction based on:
- a) RNA splicing dynamics
- b) Gene count complexity — stem cells express more genes
- c) Experimental time point labels
- d) Marker gene expression patterns

<details><summary>Answer</summary>b) CytoTRACEKernel uses the principle that less differentiated (more stem-like) cells express more genes. It estimates developmental potential from gene count complexity.</details>

---

### Q7. When would you choose CytoTRACEKernel over VelocityKernel?
- a) When velocity data is unavailable or of poor quality
- b) When you have a time-course experiment
- c) When studying non-developmental processes like T cell exhaustion
- d) When you want faster computation

<details><summary>Answer</summary>a) CytoTRACEKernel is the primary choice when RNA velocity data is not available (no spliced/unspliced counts) or when velocity quality is poor (noisy, inconsistent arrows).</details>

---

### Q8. The PseudotimeKernel accepts input from:
- a) Only diffusion pseudotime
- b) Only Monocle 3 pseudotime
- c) Any pseudotime stored in adata.obs, from any method
- d) Only CytoTRACE pseudotime

<details><summary>Answer</summary>c) PseudotimeKernel is agnostic to the pseudotime source — it accepts any numeric pseudotime column in `adata.obs`, whether from DPT, Monocle 3, Slingshot, CytoTRACE, or any custom method.</details>

---

### Q9. The RealTimeKernel connects cells across time points using:
- a) RNA velocity arrows
- b) Nearest neighbor graphs
- c) Optimal transport
- d) Gene regulatory networks

<details><summary>Answer</summary>c) RealTimeKernel uses optimal transport to infer the most likely mapping between cell populations at consecutive experimental time points, solving the problem of destructive sampling.</details>

---

### Q10. Which kernel provides the strongest directional signal?
- a) VelocityKernel — because RNA velocity captures real physical dynamics
- b) RealTimeKernel — because it uses actual measured experimental time
- c) CytoTRACEKernel — because gene counts are always reliable
- d) ConnectivityKernel — because transcriptomic similarity is fundamental

<details><summary>Answer</summary>b) RealTimeKernel uses actual experimental time, which is the most direct measurement of temporal progression. However, it requires a time-course experimental design. VelocityKernel and CytoTRACEKernel infer directionality computationally, making them useful when real time data is unavailable.</details>

---

## Part C: Kernel Combination & Analysis (Module 12)

### Q11. What is the main benefit of combining multiple kernels?
- a) It makes the computation faster
- b) It produces a more robust fate estimate by leveraging complementary data views
- c) It reduces the number of genes needed
- d) It automatically corrects batch effects

<details><summary>Answer</summary>b) Combining kernels leverages complementary data views (e.g., local splicing dynamics + global developmental potential), producing fate estimates that are more robust than any single kernel alone.</details>

---

### Q12. When two kernels strongly disagree on the fate of a cell population, you should:
- a) Always trust the VelocityKernel
- b) Discard both results
- c) Investigate biologically — the disagreement may reveal interesting biology or data quality issues
- d) Average the results and ignore the disagreement

<details><summary>Answer</summary>c) Disagreement between kernels is informative. It may reveal interesting biology (e.g., a population with ongoing transcriptional changes but stable developmental potential) or highlight data quality issues in one view. Always investigate and report.</details>

---

### Q13. In `combined = 0.6 * vk + 0.4 * ctk`, what do the weights represent?
- a) p-values for each kernel
- b) The relative contribution of each kernel's transition matrix to the combined transition matrix
- c) The number of cells used by each kernel
- d) The percentage of genes considered by each kernel

<details><summary>Answer</summary>b) The weights determine the relative contribution of each kernel's transition matrix to the combined matrix. 60% of the transition probabilities come from VelocityKernel and 40% from CytoTRACEKernel.</details>

---

### Q14. What are "driver genes" in CellRank?
- a) Genes required for cell survival
- b) Genes whose expression correlates with the fate probability toward a specific terminal state
- c) Housekeeping genes used for normalization
- d) Genes with the highest variance

<details><summary>Answer</summary>b) Driver genes are genes whose expression correlates with the absorption probability toward a specific terminal state. They are candidate regulators of lineage commitment.</details>

---

### Q15. You run CellRank with VelocityKernel and identify 3 terminal states. Then you run CytoTRACEKernel and identify 4 terminal states. What is the best next step?
- a) Only report the VelocityKernel results
- b) Only report the CytoTRACEKernel results since it found more states
- c) Report both, investigate which terminal state is unique to CytoTRACEKernel, and run a combined kernel
- d) Average the number and report 3.5 terminal states

<details><summary>Answer</summary>c) Report both results, investigate what the extra terminal state represents biologically, and run a combined kernel analysis to see which terminal states are supported by multiple data views. The additional state found by CytoTRACEKernel may be biologically meaningful or an artifact.</details>

---

## Scoring

| Score | Assessment |
|-------|------------|
| 13-15 | Excellent — ready for the final project |
| 10-12 | Good — review specific modules where you missed questions |
| 7-9 | Fair — re-read the CellRank 2 guide and redo Labs 8-9D |
| <7 | Needs review — start from Module 7 and work through all labs again |
