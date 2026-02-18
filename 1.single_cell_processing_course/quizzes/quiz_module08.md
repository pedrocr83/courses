# Quiz: Module 8 - Quality Control Metrics

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What are the three core per-cell QC metrics?

- A) Gene count, cell count, batch ID
- B) Number of genes detected, total UMI counts, mitochondrial percentage
- C) Mapping rate, duplication rate, error rate
- D) Mean expression, variance, sparsity

**Your answer:** _____

---

### Question 2
A cell with 25% mitochondrial reads likely indicates:

- A) A healthy, high-quality cell
- B) A dying or stressed cell
- C) A doublet
- D) An empty droplet

**Your answer:** _____

---

### Question 3
What does a very high gene count per cell often suggest?

- A) High-quality cell
- B) Possible doublet (two cells captured as one)
- C) Dead cell
- D) Empty droplet

**Your answer:** _____

---

### Question 4
Why is mitochondrial content used as a QC metric?

- A) Mitochondrial genes are always noise
- B) Dying cells release cytoplasmic RNA but retain mitochondrial RNA
- C) Healthy cells have no mitochondrial genes
- D) It's required by Cell Ranger

**Your answer:** _____

---

### Question 5
What is the typical mitochondrial percentage threshold for most tissues?

- A) 0-1%
- B) 5-20%
- C) 50-80%
- D) There is no standard - it's tissue-dependent

**Your answer:** _____

---

### Question 6
Per-gene QC typically filters out genes expressed in:

- A) More than 50% of cells
- B) Fewer than 3 cells
- C) All cells
- D) Only mitochondrial genes

**Your answer:** _____

---

### Question 7
What visualization is best for examining QC metric distributions?

- A) Heatmap
- B) Violin plot or histogram
- C) Network graph
- D) Phylogenetic tree

**Your answer:** _____

---

### Question 8
A scatter plot of genes vs UMIs per cell helps identify:

- A) Batch effects
- B) Cells with unusual gene/UMI ratios (potential quality issues)
- C) Cell types
- D) Differential expression

**Your answer:** _____

---

### Question 9
What does MAD stand for in the context of adaptive thresholds?

- A) Maximum Allowed Deviation
- B) Median Absolute Deviation
- C) Mean Adjusted Distance
- D) Minimum Acceptable Density

**Your answer:** _____

---

### Question 10
Why should QC thresholds be dataset-specific?

- A) To make analysis faster
- B) Different tissues and conditions have different expected values
- C) Software requires it
- D) For publication requirements

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
Explain why a dying cell might have high mitochondrial percentage but low total UMI count.

**Your answer:**

---

### Question 12
What is the difference between a fixed threshold and an adaptive (MAD-based) threshold for QC?

**Your answer:**

---

### Question 13
Why might brain tissue have a higher acceptable mitochondrial percentage than blood?

**Your answer:**

---

### Question 14
What does a gene vs UMI scatter plot show that violin plots of individual metrics might miss?

**Your answer:**

---

### Question 15
Give one reason why you might NOT filter out genes expressed in fewer than 3 cells.

**Your answer:**

---

## Answer Key

### Part A
1. B  
2. B  
3. B  
4. B  
5. D  
6. B  
7. B  
8. B  
9. B  
10. B  

### Part B (sample answers)
11. Dying cells lose cytoplasmic (nuclear-encoded) mRNA due to leakage/RNase, while mitochondrial mRNA is retained because mitochondria remain intact longer.  
12. Fixed: use a single cutoff (e.g. mt% < 10%). Adaptive: use MAD to define outliers relative to the dataset's distribution, adapting to tissue/condition.  
13. Neurons have many mitochondria (high energy demand); brain tissue naturally has higher mt% than blood cells.  
14. Joint distribution – cells with high genes AND high UMIs may be doublets; low genes + low UMIs may be empty/dying; the ratio reveals patterns.  
15. Rare cell types or rare genes of biological interest; removing them could eliminate important biology before analysis.
