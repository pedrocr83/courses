# Quiz: Module 11 - Pseudobulk Analysis

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What is pseudobulk analysis?

- A) Analyzing bulk RNA-seq as if it were single-cell
- B) Aggregating single-cell counts to sample level
- C) Removing bulk samples from analysis
- D) A type of normalization

**Your answer:** _____

---

### Question 2
Why is pseudobulk recommended for scRNA-seq DE?

- A) It's faster to compute
- B) It respects the biological replicate structure
- C) It detects more genes
- D) It requires fewer cells

**Your answer:** _____

---

### Question 3
What is "pseudoreplication" in scRNA-seq DE?

- A) Having too few replicates
- B) Treating cells as independent when they come from the same sample
- C) Duplicating samples
- D) Technical replicates

**Your answer:** _____

---

### Question 4
How are cells typically aggregated for pseudobulk?

- A) By gene
- B) By sample/donor and cell type
- C) By expression level
- D) Randomly

**Your answer:** _____

---

### Question 5
After pseudobulk aggregation, which tools can be used for DE?

- A) Only Seurat
- B) Only Scanpy
- C) Bulk RNA-seq tools like DESeq2 and edgeR
- D) Only MAST

**Your answer:** _____

---

### Question 6
Pseudobulk analysis requires that your experiment has:

- A) At least 10,000 cells
- B) Multiple biological replicates per condition
- C) Paired samples
- D) Only one cell type

**Your answer:** _____

---

### Question 7
What is the main limitation of cell-level Wilcoxon tests?

- A) Too slow
- B) Inflate false positives due to pseudoreplication
- C) Cannot detect upregulated genes
- D) Require normalization

**Your answer:** _____

---

### Question 8
When might cell-level tests be appropriate?

- A) For publication-quality DE
- B) For exploratory marker gene identification
- C) When you have many biological replicates
- D) Never

**Your answer:** _____

---

### Question 9
In pseudobulk, what happens if a cell type has few cells in a sample?

- A) Those cells are removed
- B) Counts may be unreliable; consider filtering
- C) Results are automatically corrected
- D) The sample is duplicated

**Your answer:** _____

---

### Question 10
Pseudobulk DE typically identifies:

- A) More DEGs than cell-level tests
- B) Fewer but more reliable DEGs
- C) Exactly the same DEGs
- D) Only marker genes

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
Explain in 1–2 sentences why treating each cell as an independent replicate inflates false positives in DE.

**Your answer:**

---

### Question 12
What aggregation strategy (sum, mean, or weighted) is most common for pseudobulk, and why?

**Your answer:**

---

### Question 13
If you have 3 donors per condition and 2 cell types, how many pseudobulk "samples" do you have for a cell-type-specific DE analysis?

**Your answer:**

---

### Question 14
When would you NOT use pseudobulk for scRNA-seq DE?

**Your answer:**

---

### Question 15
What is the minimum recommended number of cells per sample/cell-type for reliable pseudobulk?

**Your answer:**

---

## Answer Key

### Part A
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

### Part B (sample answers)
11. Cells from the same donor are correlated; treating them as independent inflates the effective sample size, leading to artificially small p-values and false positives.  
12. Sum (aggregating counts) is most common; preserves count structure for DESeq2/edgeR; maintains interpretability.  
13. 3 donors × 2 conditions × 2 cell types = 12 pseudobulk samples (one per donor–condition–cell type).  
14. When you have only one biological replicate per condition; for exploratory marker discovery; when sample is too small.  
15. Typically 10–50+ cells per sample/cell-type; some recommend 100+ for stable estimates.
