# Quiz: Module 2 - Principal Component Analysis

**10 Questions | Passing: 70%**

---

## Question 1
What does PCA do to high-dimensional data?

- A) Increases the number of dimensions
- B) Finds linear combinations that maximize variance
- C) Clusters cells automatically
- D) Normalizes gene expression

---

## Question 2
The first principal component (PC1) captures:

- A) The least variance in the data
- B) The most variance in the data
- C) Exactly 50% of variance
- D) Only technical variation

---

## Question 3
Why do we select "highly variable genes" before PCA?

- A) To speed up computation only
- B) To focus on biologically informative genes and reduce noise
- C) PCA requires exactly 2000 genes
- D) To remove all housekeeping genes

---

## Question 4
The "elbow plot" helps you decide:

- A) Which genes to remove
- B) How many principal components to use
- C) The clustering resolution
- D) Cell type annotations

---

## Question 5
What do PC loadings tell you?

- A) How many cells are in each cluster
- B) Which genes contribute most to each PC
- C) The p-value of each gene
- D) Cell type proportions

---

## Question 6
If PC1 is highly correlated with total UMI counts, this suggests:

- A) Perfect normalization
- B) Incomplete normalization / technical variation
- C) The data has no biological signal
- D) Too many PCs were selected

---

## Question 7
Typical number of PCs used for downstream analysis in scRNA-seq:

- A) 1-2
- B) 10-50
- C) 500-1000
- D) All PCs

---

## Question 8
PCA assumes that the data structure is:

- A) Circular
- B) Linear
- C) Branching
- D) Random

---

## Question 9
What happens if you use too few PCs?

- A) Lose biological variation, miss rare cell types
- B) Include too much noise
- C) Faster but identical results
- D) Clustering will fail completely

---

## Question 10
Before running PCA on scRNA-seq data, you should:

- A) Run clustering first
- B) Log-normalize and select HVGs
- C) Convert to protein data
- D) Remove all mitochondrial genes

---

## Answer Key

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

