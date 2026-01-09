# Quiz: Module 11 - Pseudobulk Analysis

**10 Questions | Passing: 70%**

---

## Question 1
What is pseudobulk analysis?

- A) Analyzing bulk RNA-seq as if it were single-cell
- B) Aggregating single-cell counts to sample level
- C) Removing bulk samples from analysis
- D) A type of normalization

---

## Question 2
Why is pseudobulk recommended for scRNA-seq DE?

- A) It's faster to compute
- B) It respects the biological replicate structure
- C) It detects more genes
- D) It requires fewer cells

---

## Question 3
What is "pseudoreplication" in scRNA-seq DE?

- A) Having too few replicates
- B) Treating cells as independent when they come from the same sample
- C) Duplicating samples
- D) Technical replicates

---

## Question 4
How are cells typically aggregated for pseudobulk?

- A) By gene
- B) By sample/donor and cell type
- C) By expression level
- D) Randomly

---

## Question 5
After pseudobulk aggregation, which tools can be used for DE?

- A) Only Seurat
- B) Only Scanpy
- C) Bulk RNA-seq tools like DESeq2 and edgeR
- D) Only MAST

---

## Question 6
Pseudobulk analysis requires that your experiment has:

- A) At least 10,000 cells
- B) Multiple biological replicates per condition
- C) Paired samples
- D) Only one cell type

---

## Question 7
What is the main limitation of cell-level Wilcoxon tests?

- A) Too slow
- B) Inflate false positives due to pseudoreplication
- C) Cannot detect upregulated genes
- D) Require normalization

---

## Question 8
When might cell-level tests be appropriate?

- A) For publication-quality DE
- B) For exploratory marker gene identification
- C) When you have many biological replicates
- D) Never

---

## Question 9
In pseudobulk, what happens if a cell type has few cells in a sample?

- A) Those cells are removed
- B) Counts may be unreliable; consider filtering
- C) Results are automatically corrected
- D) The sample is duplicated

---

## Question 10
Pseudobulk DE typically identifies:

- A) More DEGs than cell-level tests
- B) Fewer but more reliable DEGs
- C) Exactly the same DEGs
- D) Only marker genes

---

## Answer Key

1. B
2. B
3. B
4. B
5. C
6. B
7. B
8. B
9. B
10. B

