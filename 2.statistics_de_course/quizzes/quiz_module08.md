# Quiz: Module 8 - Bulk RNA-seq DE with DESeq2

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What statistical model does DESeq2 use?

- A) Linear regression
- B) Negative binomial generalized linear model
- C) Poisson regression
- D) Logistic regression

**Your answer:** _____

---

### Question 2
What is the purpose of DESeq2's size factors?

- A) Determine statistical significance
- B) Account for differences in sequencing depth
- C) Calculate fold changes
- D) Remove batch effects

**Your answer:** _____

---

### Question 3
What does DESeq2's shrinkage of fold changes accomplish?

- A) Increases all fold changes
- B) Reduces noisy estimates for low-count genes
- C) Removes all non-significant genes
- D) Normalizes the data

**Your answer:** _____

---

### Question 4
In DESeq2, what does the `baseMean` column represent?

- A) p-value
- B) Average normalized count across all samples
- C) Fold change
- D) Dispersion estimate

**Your answer:** _____

---

### Question 5
What does `padj` in DESeq2 results represent?

- A) Raw p-value
- B) FDR-adjusted p-value (Benjamini-Hochberg)
- C) Bonferroni-corrected p-value
- D) Effect size

**Your answer:** _____

---

### Question 6
Which function runs the full DESeq2 analysis pipeline?

- A) `DESeqDataSetFromMatrix()`
- B) `DESeq()`
- C) `results()`
- D) `estimateSizeFactors()`

**Your answer:** _____

---

### Question 7
What is the recommended minimum number of biological replicates for DESeq2?

- A) 1
- B) 2
- C) 3 or more
- D) 10

**Your answer:** _____

---

### Question 8
When should you include batch in the DESeq2 design formula?

- A) Never
- B) Always, regardless of design
- C) When samples were processed in different batches
- D) Only for single-cell data

**Your answer:** _____

---

### Question 9
What does a positive log2FoldChange indicate?

- A) Gene is downregulated in the test condition
- B) Gene is upregulated in the test condition
- C) Gene has low expression
- D) Result is not significant

**Your answer:** _____

---

### Question 10
Before running DESeq2, you should filter out genes with:

- A) High expression
- B) Negative values
- C) Very low counts across samples
- D) Large fold changes

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
What is the design formula in DESeq2, and what does `~ condition + batch` mean?

**Your answer:**

---

### Question 12
Why does DESeq2 shrink dispersion estimates? What problem does this address?

**Your answer:**

---

### Question 13
If your design is `~ condition` and you have 3 conditions (A, B, C), how does DESeq2 handle the comparison?

**Your answer:**

---

### Question 14
Why filter low-count genes before DESeq2 analysis?

**Your answer:**

---

### Question 15
What is the difference between `pvalue` and `padj` in DESeq2 results, and when would you use each?

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
7. C  
8. C  
9. B  
10. C  

### Part B (sample answers)
11. Design formula specifies the model. `~ condition + batch` means expression is modeled by condition and batch; batch is a nuisance covariate to adjust for.  
12. Shrinkage borrows information across genes to get better dispersion estimates for genes with few replicates; avoids unstable estimates and inflated/deflated p-values.  
13. DESeq2 estimates coefficients for each condition; you specify contrasts (e.g. B vs A, C vs A) in `results()` to get comparisons.  
14. Low-count genes have highly variable estimates and low power; they add to multiple testing burden and often produce unreliable results.  
15. `pvalue`: raw p-value. `padj`: FDR-adjusted; use for calling DEGs to control false discoveries. Use `pvalue` for diagnostic plots.
