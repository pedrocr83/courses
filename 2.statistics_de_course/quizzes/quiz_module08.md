# Quiz: Module 8 - Bulk RNA-seq DE with DESeq2

**10 Questions | Passing: 70%**

---

## Question 1
What statistical model does DESeq2 use?

- A) Linear regression
- B) Negative binomial generalized linear model
- C) Poisson regression
- D) Logistic regression

---

## Question 2
What is the purpose of DESeq2's size factors?

- A) Determine statistical significance
- B) Account for differences in sequencing depth
- C) Calculate fold changes
- D) Remove batch effects

---

## Question 3
What does DESeq2's shrinkage of fold changes accomplish?

- A) Increases all fold changes
- B) Reduces noisy estimates for low-count genes
- C) Removes all non-significant genes
- D) Normalizes the data

---

## Question 4
In DESeq2, what does the `baseMean` column represent?

- A) p-value
- B) Average normalized count across all samples
- C) Fold change
- D) Dispersion estimate

---

## Question 5
What does `padj` in DESeq2 results represent?

- A) Raw p-value
- B) FDR-adjusted p-value (Benjamini-Hochberg)
- C) Bonferroni-corrected p-value
- D) Effect size

---

## Question 6
Which function runs the full DESeq2 analysis pipeline?

- A) `DESeqDataSetFromMatrix()`
- B) `DESeq()`
- C) `results()`
- D) `estimateSizeFactors()`

---

## Question 7
What is the recommended minimum number of biological replicates for DESeq2?

- A) 1
- B) 2
- C) 3 or more
- D) 10

---

## Question 8
When should you include batch in the DESeq2 design formula?

- A) Never
- B) Always, regardless of design
- C) When samples were processed in different batches
- D) Only for single-cell data

---

## Question 9
What does a positive log2FoldChange indicate?

- A) Gene is downregulated in the test condition
- B) Gene is upregulated in the test condition
- C) Gene has low expression
- D) Result is not significant

---

## Question 10
Before running DESeq2, you should filter out genes with:

- A) High expression
- B) Negative values
- C) Very low counts across samples
- D) Large fold changes

---

## Answer Key

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

