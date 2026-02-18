# Quiz: Module 9 - CellRank

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
CellRank models cell state transitions as:

- A) Linear regression
- B) Markov chain / transition probabilities
- C) Neural network
- D) Clustering

**Your answer:** _____

---

### Question 2
CellRank can combine:

- A) Only RNA velocity
- B) RNA velocity with transcriptomic similarity
- C) Only gene expression
- D) Protein data only

**Your answer:** _____

---

### Question 3
Terminal states in CellRank represent:

- A) Dead cells
- B) Final differentiated cell fates
- C) Starting cells
- D) Technical artifacts

**Your answer:** _____

---

### Question 4
Fate probabilities in CellRank tell you:

- A) Gene expression levels
- B) The likelihood each cell reaches each terminal state
- C) Cluster membership
- D) Batch effects

**Your answer:** _____

---

### Question 5
Absorption probabilities are calculated using:

- A) Simple averaging
- B) Solving the Markov chain fundamental matrix
- C) Random sampling
- D) Gradient descent

**Your answer:** _____

---

### Question 6
CellRank's "macrostates" are:

- A) Individual cells
- B) Metastable cell populations (grouped states)
- C) Gene modules
- D) Technical batches

**Your answer:** _____

---

### Question 7
The transition matrix in CellRank has dimensions:

- A) genes × cells
- B) cells × cells
- C) genes × genes
- D) clusters × clusters

**Your answer:** _____

---

### Question 8
CellRank is particularly useful for:

- A) Simple two-state systems
- B) Multi-lineage differentiation with branch points
- C) Data without any trajectory structure
- D) Bulk RNA-seq only

**Your answer:** _____

---

### Question 9
"Driver genes" in CellRank are:

- A) Housekeeping genes
- B) Genes highly correlated with fate probabilities
- C) The most expressed genes
- D) Ribosomal genes

**Your answer:** _____

---

### Question 10
Compared to simple pseudotime, CellRank provides:

- A) Only ordering
- B) Directional fate probabilities and driver genes
- C) Faster computation
- D) Simpler interpretation

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
What is a transition matrix in the context of CellRank, and what does each entry represent?

**Your answer:**

---

### Question 12
How does CellRank use RNA velocity to build the transition matrix?

**Your answer:**

---

### Question 13
What is the difference between a "terminal state" and a "macrostate"?

**Your answer:**

---

### Question 14
Why are driver genes useful for biological interpretation?

**Your answer:**

---

### Question 15
When might you combine VelocityKernel with ConnectivityKernel in CellRank 2, and why?

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
11. cells × cells matrix; entry (i,j) = probability of transitioning from cell i to cell j (or cell j's neighborhood).  
12. Velocity gives direction; higher velocity toward j increases transition probability from i to j.  
13. Terminal: absorbing state (differentiation endpoint). Macrostate: metastable group of cells (GPCCA clusters).  
14. They are candidate regulators of fate commitment; expression correlates with probability of reaching a given fate.  
15. Velocity is directional but noisy; connectivity smooths; combine (e.g. 0.8*vk + 0.2*ck) for robustness.
