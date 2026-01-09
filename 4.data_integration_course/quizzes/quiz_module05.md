# Quiz: Module 5 - Harmony Integration

**10 Questions | Passing: 70%**

---

## Question 1
Harmony operates on which data representation?

- A) Raw counts
- B) PCA embedding
- C) Gene expression matrix
- D) UMAP coordinates

---

## Question 2
Harmony uses what type of clustering approach?

- A) Hard clustering (one cluster per cell)
- B) Soft clustering (cells belong to multiple clusters with probabilities)
- C) No clustering
- D) Hierarchical clustering

---

## Question 3
The main advantage of Harmony is:

- A) Highest accuracy on all datasets
- B) Speed and scalability to large datasets
- C) Only works with Seurat
- D) Requires no parameters

---

## Question 4
Harmony iteratively:

- A) Removes cells
- B) Adjusts embeddings to mix batches within soft clusters
- C) Changes gene expression
- D) Selects variable genes

---

## Question 5
After Harmony integration, you should:

- A) Re-run PCA
- B) Use the corrected embedding for clustering and UMAP
- C) Delete the original data
- D) Skip all downstream analysis

---

## Question 6
The "theta" parameter in Harmony controls:

- A) Number of iterations
- B) Diversity penalty for batch mixing
- C) Learning rate
- D) Number of clusters

---

## Question 7
Harmony preserves biological variation by:

- A) Removing all technical variation
- B) Only correcting within similar cell types (soft clusters)
- C) Averaging all cells together
- D) Using raw counts

---

## Question 8
Which function runs Harmony in Seurat?

- A) `FindClusters()`
- B) `RunHarmony()`
- C) `IntegrateData()`
- D) `NormalizeData()`

---

## Question 9
Which function runs Harmony in Scanpy?

- A) `sc.pp.harmony()`
- B) `sc.external.pp.harmony_integrate()`
- C) `sc.tl.harmony()`
- D) `sce.harmonize()`

---

## Question 10
Harmony is best suited for:

- A) Only 2 batches
- B) Multiple batches with shared cell types
- C) Data without any batch effects
- D) Proteomics data only

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
9. B
10. B

