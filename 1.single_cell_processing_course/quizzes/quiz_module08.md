# Quiz: Module 8 - Quality Control Metrics

**10 Questions | Passing: 70%**

---

## Question 1
What are the three core per-cell QC metrics?

- A) Gene count, cell count, batch ID
- B) Number of genes detected, total UMI counts, mitochondrial percentage
- C) Mapping rate, duplication rate, error rate
- D) Mean expression, variance, sparsity

---

## Question 2
A cell with 25% mitochondrial reads likely indicates:

- A) A healthy, high-quality cell
- B) A dying or stressed cell
- C) A doublet
- D) An empty droplet

---

## Question 3
What does a very high gene count per cell often suggest?

- A) High-quality cell
- B) Possible doublet (two cells captured as one)
- C) Dead cell
- D) Empty droplet

---

## Question 4
Why is mitochondrial content used as a QC metric?

- A) Mitochondrial genes are always noise
- B) Dying cells release cytoplasmic RNA but retain mitochondrial RNA
- C) Healthy cells have no mitochondrial genes
- D) It's required by Cell Ranger

---

## Question 5
What is the typical mitochondrial percentage threshold for most tissues?

- A) 0-1%
- B) 5-20%
- C) 50-80%
- D) There is no standard - it's tissue-dependent

---

## Question 6
Per-gene QC typically filters out genes expressed in:

- A) More than 50% of cells
- B) Fewer than 3 cells
- C) All cells
- D) Only mitochondrial genes

---

## Question 7
What visualization is best for examining QC metric distributions?

- A) Heatmap
- B) Violin plot or histogram
- C) Network graph
- D) Phylogenetic tree

---

## Question 8
A scatter plot of genes vs UMIs per cell helps identify:

- A) Batch effects
- B) Cells with unusual gene/UMI ratios (potential quality issues)
- C) Cell types
- D) Differential expression

---

## Question 9
What does MAD stand for in the context of adaptive thresholds?

- A) Maximum Allowed Deviation
- B) Median Absolute Deviation
- C) Mean Adjusted Distance
- D) Minimum Acceptable Density

---

## Question 10
Why should QC thresholds be dataset-specific?

- A) To make analysis faster
- B) Different tissues and conditions have different expected values
- C) Software requires it
- D) For publication requirements

---

## Answer Key

1. B
2. B
3. B
4. B
5. D (but 5-20% is a common starting point)
6. B
7. B
8. B
9. B
10. B

