# Quiz: Module 8 - Evaluating Integration Quality

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
The two main aspects of integration quality are:

- A) Speed and memory
- B) Batch mixing and biological conservation
- C) Number of cells and genes
- D) PCA and UMAP

**Your answer:** _____

---

### Question 2
iLISI measures:

- A) Biological conservation
- B) Batch mixing quality
- C) Gene expression variance
- D) Cluster stability

**Your answer:** _____

---

### Question 3
cLISI measures:

- A) Batch mixing
- B) Cell type separation (biological conservation)
- C) Sequencing depth
- D) Integration speed

**Your answer:** _____

---

### Question 4
The ideal integration result has:

- A) High iLISI, High cLISI
- B) High iLISI, Low cLISI
- C) Low iLISI, High cLISI
- D) Low iLISI, Low cLISI

**Your answer:** _____

---

### Question 5
ARI (Adjusted Rand Index) compares:

- A) Gene expression levels
- B) Cluster assignments before and after integration
- C) Batch sizes
- D) Cell sizes

**Your answer:** _____

---

### Question 6
scIB is:

- A) A clustering method
- B) A benchmarking framework for integration methods
- C) A new sequencing technology
- D) A programming language

**Your answer:** _____

---

### Question 7
If integration makes cell types indistinguishable, this is:

- A) Perfect integration
- B) Over-correction
- C) Under-correction
- D) Normal behavior

**Your answer:** _____

---

### Question 8
Graph connectivity measures whether:

- A) The computer has internet
- B) Same cell type from different batches connects in the neighbor graph
- C) Genes are connected to proteins
- D) Clusters are the right size

**Your answer:** _____

---

### Question 9
A good strategy for evaluation is:

- A) Visual inspection only
- B) Multiple metrics covering both mixing and conservation
- C) Using a single metric
- D) Skipping evaluation

**Your answer:** _____

---

### Question 10
If batches remain separated after integration, you should:

- A) Declare success
- B) Try a different method or adjust parameters
- C) Remove one batch
- D) Add more batches

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
Why must integration evaluation consider both batch mixing AND biological conservation?

**Your answer:**

---

### Question 12
What does graph connectivity tell you, and why does it matter for downstream analysis?

**Your answer:**

---

### Question 13
How would you use known marker genes to validate integration quality?

**Your answer:**

---

### Question 14
What is the trade-off between batch mixing and biological conservation?

**Your answer:**

---

### Question 15
If two integration methods give conflicting metrics (e.g., one better iLISI, other better cLISI), how would you decide?

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
11. Mixing alone: over-correction can merge distinct cell types. Conservation alone: under-correction leaves batch structure. Both needed.  
12. Whether cells of same type from different batches are neighbors; affects clustering and trajectory inference.  
13. Check that known cell-type markers still separate cell types; expression patterns preserved.  
14. Aggressive mixing may erase real biological differences; weak mixing leaves batch structure. Balance is key.  
15. Consider biological goals; use multiple metrics; validate with known biology; report both and discuss.
