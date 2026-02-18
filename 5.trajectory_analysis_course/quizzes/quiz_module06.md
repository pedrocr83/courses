# Quiz: Module 6 - RNA Velocity

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
RNA velocity uses the ratio of:

- A) Expressed vs non-expressed genes
- B) Spliced vs unspliced RNA
- C) mRNA vs protein
- D) Nuclear vs cytoplasmic RNA

**Your answer:** _____

---

### Question 2
Unspliced RNA contains:

- A) Only exons
- B) Introns (not yet spliced out)
- C) Only UTRs
- D) Protein sequences

**Your answer:** _____

---

### Question 3
RNA velocity predicts:

- A) Past cell states
- B) Future cell states / direction of change
- C) Cell cycle phase only
- D) Mutation rates

**Your answer:** _____

---

### Question 4
Velocity arrows on a UMAP point toward:

- A) Higher gene expression
- B) Predicted future cell state
- C) The nearest neighbor
- D) Technical artifacts

**Your answer:** _____

---

### Question 5
To compute RNA velocity, you need:

- A) Only a count matrix
- B) Separate spliced and unspliced count matrices
- C) Protein abundance data
- D) ChIP-seq data

**Your answer:** _____

---

### Question 6
scVelo's "stochastic" model assumes:

- A) No noise in the data
- B) Constant transcription/splicing rates
- C) Only one gene is expressed
- D) Circular trajectories

**Your answer:** _____

---

### Question 7
scVelo's "dynamical" model:

- A) Is faster than stochastic
- B) Learns gene-specific kinetic parameters
- C) Requires labeled data
- D) Only works on small datasets

**Your answer:** _____

---

### Question 8
Chaotic/noisy velocity arrows suggest:

- A) Perfect differentiation trajectory
- B) Steady state cells or poor velocity estimation
- C) Too much data
- D) Correct root selection

**Your answer:** _____

---

### Question 9
Velocity data is typically generated using:

- A) Standard Cell Ranger
- B) velocyto or kallisto | bustools (lamanno workflow)
- C) FASTQC
- D) Manual counting

**Your answer:** _____

---

### Question 10
A limitation of RNA velocity is:

- A) Only works for 10x data
- B) Assumes specific splicing kinetics that may not hold for all genes
- C) Requires time-series data
- D) Cannot visualize results

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
In one sentence, explain how spliced/unspliced counts inform the "direction" of expression change.

**Your answer:**

---

### Question 12
What is the difference between the stochastic and dynamical models in scVelo?

**Your answer:**

---

### Question 13
Why might velocity arrows be unreliable in steady-state cell populations?

**Your answer:**

---

### Question 14
What preprocessing or quantification is required to run RNA velocity (different from standard 3' scRNA-seq)?

**Your answer:**

---

### Question 15
Name one biological scenario where RNA velocity is particularly useful vs one where it may fail.

**Your answer:**

---

## Answer Key

### Part A
1. B  
2. B  
3. B  
4. B  
5. B  
6. B  
7. B  
8. B  
9. B  
10. B  

### Part B (sample answers)
11. More unspliced than spliced → gene is being induced; more spliced → gene is being repressed; ratio gives direction.  
12. Stochastic: constant rates, faster. Dynamical: learns gene-specific kinetics from data, more accurate but slower.  
13. In steady state, spliced/unspliced are balanced; no net direction; arrows become noisy or meaningless.  
14. Need spliced and unspliced counts; typically 3' data with intronic reads (velocyto, kallisto|bustools lamanno).  
15. Useful: differentiation, activation. May fail: steady-state cells, genes with nonstandard splicing, poor quality data.
