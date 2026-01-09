# Quiz: Module 2 - Diagnosing Batch Effects

**10 Questions | Passing: 70%**

---

## Question 1
What is a batch effect?

- A) Biological differences between cell types
- B) Technical variation due to experimental processing
- C) Random noise
- D) Sequencing errors only

---

## Question 2
Which visualization best reveals batch effects?

- A) Violin plot of gene expression
- B) UMAP colored by batch/sample
- C) Heatmap of top genes
- D) Bar chart of cell counts

---

## Question 3
LISI (Local Inverse Simpson Index) measures:

- A) Gene expression levels
- B) Diversity of batch labels in local neighborhoods
- C) Sequencing quality
- D) Number of clusters

---

## Question 4
High iLISI (integration LISI) indicates:

- A) Poor batch mixing
- B) Good batch mixing
- C) Strong biological signal
- D) Low cell quality

---

## Question 5
cLISI (cell-type LISI) should be:

- A) High (mixed cell types locally)
- B) Low (same cell type locally) - conservation
- C) Equal to iLISI
- D) Zero

---

## Question 6
kBET tests whether:

- A) Genes are differentially expressed
- B) Local neighborhoods have expected batch proportions
- C) Cells are properly normalized
- D) Clustering is optimal

---

## Question 7
When is batch correction NOT needed?

- A) Never - always correct
- B) When batches are well-mixed and cell types separate
- C) When you have multiple samples
- D) For large datasets only

---

## Question 8
Over-correction of batch effects can cause:

- A) Better visualization
- B) Loss of true biological differences
- C) Faster computation
- D) More accurate annotation

---

## Question 9
Confounding occurs when:

- A) Batch and biological condition are not separable
- B) Too many cells are sequenced
- C) The wrong clustering is used
- D) Genes are lowly expressed

---

## Question 10
Before integrating, you should always:

- A) Remove all batch information
- B) QC and process each batch separately first
- C) Combine all raw data
- D) Annotate cell types

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

