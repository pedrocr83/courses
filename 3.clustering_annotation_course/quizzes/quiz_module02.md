# Quiz: Module 2 - Principal Component Analysis

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What does PCA do to high-dimensional data?

- A) Increases the number of dimensions
- B) Finds linear combinations that maximize variance
- C) Clusters cells automatically
- D) Normalizes gene expression

**Your answer:** _____

---

### Question 2
The first principal component (PC1) captures:

- A) The least variance in the data
- B) The most variance in the data
- C) Exactly 50% of variance
- D) Only technical variation

**Your answer:** _____

---

### Question 3
Why do we select "highly variable genes" before PCA?

- A) To speed up computation only
- B) To focus on biologically informative genes and reduce noise
- C) PCA requires exactly 2000 genes
- D) To remove all housekeeping genes

**Your answer:** _____

---

### Question 4
The "elbow plot" helps you decide:

- A) Which genes to remove
- B) How many principal components to use
- C) The clustering resolution
- D) Cell type annotations

**Your answer:** _____

---

### Question 5
What do PC loadings tell you?

- A) How many cells are in each cluster
- B) Which genes contribute most to each PC
- C) The p-value of each gene
- D) Cell type proportions

**Your answer:** _____

---

### Question 6
If PC1 is highly correlated with total UMI counts, this suggests:

- A) Perfect normalization
- B) Incomplete normalization / technical variation
- C) The data has no biological signal
- D) Too many PCs were selected

**Your answer:** _____

---

### Question 7
Typical number of PCs used for downstream analysis in scRNA-seq:

- A) 1-2
- B) 10-50
- C) 500-1000
- D) All PCs

**Your answer:** _____

---

### Question 8
PCA assumes that the data structure is:

- A) Circular
- B) Linear
- C) Branching
- D) Random

**Your answer:** _____

---

### Question 9
What happens if you use too few PCs?

- A) Lose biological variation, miss rare cell types
- B) Include too much noise
- C) Faster but identical results
- D) Clustering will fail completely

**Your answer:** _____

---

### Question 10
Before running PCA on scRNA-seq data, you should:

- A) Run clustering first
- B) Log-normalize and select HVGs
- C) Convert to protein data
- D) Remove all mitochondrial genes

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
Why is PCA performed on highly variable genes rather than all genes?

**Your answer:**

---

### Question 12
What does "variance explained" mean for a principal component, and how is it calculated?

**Your answer:**

---

### Question 13
If PC1 explains 15% of variance and PC2 explains 8%, what does that suggest about the data structure?

**Your answer:**

---

### Question 14
How can you check if a PC captures technical (e.g., batch) vs biological variation?

**Your answer:**

---

### Question 15
What is the relationship between PCA and UMAP/t-SNE? Why is PCA typically run first?

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
9. A  
10. B  

### Part B (sample answers)
11. HVGs capture cell-to-cell differences; constant/low-variance genes add noise; reduces dimensionality and improves signal.  
12. Fraction of total variance in the data captured by that PC; sum of squared loadings / total variance.  
13. Variance is spread across many dimensions; no single dominant structure; may need many PCs.  
14. Correlate PC scores with batch, total counts, mt%; if strong correlation, it's likely technical.  
15. UMAP/t-SNE run on PCA embedding (e.g. first 30-50 PCs) for computational efficiency and noise reduction; PCA denoises before nonlinear viz.
