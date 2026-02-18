# Quiz: CellRank 2 — Multiview Fate Mapping

**Covers Modules 7-12 | 15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: CellRank 2 Framework (5 questions)

### Question 1
What is a "kernel" in CellRank 2?

- A) A machine learning algorithm for classification
- B) An object that computes a cell-cell transition matrix from a specific data view
- C) A gene expression clustering method
- D) A Python package for deep learning

**Your answer:** _____

---

### Question 2
What does GPCCA stand for, and what does it do?

- A) General Purpose Cell Cluster Analysis — clusters cells
- B) Generalized Perron Cluster Cluster Analysis — identifies metastable macrostates from the transition matrix
- C) Gene-based Principal Component Correlation Analysis — reduces dimensionality
- D) Graph-Partitioned Cell Classification Algorithm — classifies cell types

**Your answer:** _____

---

### Question 3
In CellRank 2, what are "fate probabilities"?

- A) The probability that a gene is expressed in a cell
- B) The absorption probability of each cell reaching each terminal state
- C) The likelihood that a cell type exists in the sample
- D) The p-value from differential expression testing

**Your answer:** _____

---

### Question 4
When combining VelocityKernel with ConnectivityKernel, what weighting is typically recommended?

- A) Equal weights (0.5 + 0.5)
- B) Higher weight to VelocityKernel (e.g. 0.8 velocity + 0.2 connectivity)
- C) Only ConnectivityKernel
- D) Multiply them: vk * ck

**Your answer:** _____

---

### Question 5
How do you determine the appropriate number of macrostates in GPCCA?

- A) Always use n_states=2
- B) Use the number of clusters from Leiden
- C) Examine the eigenvalue spectrum gap from Schur decomposition
- D) Count the number of genes

**Your answer:** _____

---

## Part B: Individual Kernels (5 questions)

### Question 6
CytoTRACEKernel infers developmental direction based on:

- A) RNA splicing dynamics
- B) Gene count complexity — stem cells express more genes
- C) Experimental time point labels
- D) Marker gene expression patterns

**Your answer:** _____

---

### Question 7
When would you choose CytoTRACEKernel over VelocityKernel?

- A) When velocity data is unavailable or of poor quality
- B) When you have a time-course experiment
- C) When studying non-developmental processes only
- D) When you want faster computation

**Your answer:** _____

---

### Question 8
PseudotimeKernel accepts input from:

- A) Only diffusion pseudotime
- B) Only Monocle 3 pseudotime
- C) Any pseudotime stored in adata.obs, from any method
- D) Only CytoTRACE pseudotime

**Your answer:** _____

---

### Question 9
RealTimeKernel connects cells across time points using:

- A) RNA velocity arrows
- B) Nearest neighbor graphs
- C) Optimal transport
- D) Gene regulatory networks

**Your answer:** _____

---

### Question 10
Which kernel provides the strongest directional signal when available?

- A) VelocityKernel
- B) RealTimeKernel — uses actual experimental time
- C) CytoTRACEKernel
- D) ConnectivityKernel

**Your answer:** _____

---

## Part C: Kernel Combination & Analysis (5 questions)

### Question 11
What is the main benefit of combining multiple kernels?

- A) Faster computation
- B) More robust fate estimates by leveraging complementary data views
- C) Fewer genes needed
- D) Automatic batch correction

**Your answer:** _____

---

### Question 12
When two kernels strongly disagree on the fate of a cell population, you should:

- A) Always trust VelocityKernel
- B) Discard both results
- C) Investigate biologically — disagreement may reveal biology or data quality issues
- D) Average and ignore

**Your answer:** _____

---

### Question 13
In `combined = 0.6 * vk + 0.4 * ctk`, what do the weights represent?

- A) p-values for each kernel
- B) Relative contribution of each kernel's transition matrix to the combined matrix
- C) Number of cells used
- D) Percentage of genes considered

**Your answer:** _____

---

### Question 14
What are "driver genes" in CellRank?

- A) Genes required for survival
- B) Genes whose expression correlates with fate probability toward a terminal state
- C) Housekeeping genes
- D) Highest variance genes

**Your answer:** _____

---

### Question 15
VelocityKernel finds 3 terminal states; CytoTRACEKernel finds 4. Best next step?

- A) Report only VelocityKernel
- B) Report only CytoTRACEKernel
- C) Report both, investigate the extra state, run combined kernel
- D) Average to 3.5 states

**Your answer:** _____

---

## Answer Key

### Part A
1. B  
2. B  
3. B  
4. B  
5. C  

### Part B
6. B  
7. A  
8. C  
9. C  
10. B  

### Part C
11. B  
12. C  
13. B  
14. B  
15. C  
