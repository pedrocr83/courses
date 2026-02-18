# Quiz: Module 5 - Harmony Integration

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
Harmony operates on which data representation?

- A) Raw counts
- B) PCA embedding
- C) Gene expression matrix
- D) UMAP coordinates

**Your answer:** _____

---

### Question 2
Harmony uses what type of clustering approach?

- A) Hard clustering (one cluster per cell)
- B) Soft clustering (cells belong to multiple clusters with probabilities)
- C) No clustering
- D) Hierarchical clustering

**Your answer:** _____

---

### Question 3
The main advantage of Harmony is:

- A) Highest accuracy on all datasets
- B) Speed and scalability to large datasets
- C) Only works with Seurat
- D) Requires no parameters

**Your answer:** _____

---

### Question 4
Harmony iteratively:

- A) Removes cells
- B) Adjusts embeddings to mix batches within soft clusters
- C) Changes gene expression
- D) Selects variable genes

**Your answer:** _____

---

### Question 5
After Harmony integration, you should:

- A) Re-run PCA
- B) Use the corrected embedding for clustering and UMAP
- C) Delete the original data
- D) Skip all downstream analysis

**Your answer:** _____

---

### Question 6
The "theta" parameter in Harmony controls:

- A) Number of iterations
- B) Diversity penalty for batch mixing
- C) Learning rate
- D) Number of clusters

**Your answer:** _____

---

### Question 7
Harmony preserves biological variation by:

- A) Removing all technical variation
- B) Only correcting within similar cell types (soft clusters)
- C) Averaging all cells together
- D) Using raw counts

**Your answer:** _____

---

### Question 8
Which function runs Harmony in Seurat?

- A) `FindClusters()`
- B) `RunHarmony()`
- C) `IntegrateData()`
- D) `NormalizeData()`

**Your answer:** _____

---

### Question 9
Which function runs Harmony in Scanpy?

- A) `sc.pp.harmony()`
- B) `sc.external.pp.harmony_integrate()`
- C) `sc.tl.harmony()`
- D) `sce.harmonize()`

**Your answer:** _____

---

### Question 10
Harmony is best suited for:

- A) Only 2 batches
- B) Multiple batches with shared cell types
- C) Data without any batch effects
- D) Proteomics data only

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
Why does Harmony use soft clustering rather than hard clustering?

**Your answer:**

---

### Question 12
What input does Harmony require (e.g., from Seurat or Scanpy workflow)?

**Your answer:**

---

### Question 13
How does Harmony avoid over-correction compared to methods that aggressively remove batch?

**Your answer:**

---

### Question 14
When might you increase vs decrease the theta parameter in Harmony?

**Your answer:**

---

### Question 15
What is the main limitation of Harmony compared to deep-learning methods like scVI?

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
11. Soft clustering allows cells to belong to multiple clusters; better for mixed populations and avoids arbitrary boundaries.  
12. PCA embedding and batch labels; sometimes cell type labels if using supervised mode.  
13. It corrects only within clusters of similar cells; distinct cell types stay separate.  
14. Increase theta: stronger correction when batch effects are severe. Decrease: when batch is mild or risk of over-correction.  
15. Harmony is linear (operates on PCA); scVI captures nonlinear structure and may handle complex batch effects better.
