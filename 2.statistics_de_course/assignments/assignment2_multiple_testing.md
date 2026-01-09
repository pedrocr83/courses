# Assignment 2: Multiple Testing Simulation

**Module:** 6 - Multiple Testing Correction  
**Due:** End of Week 2  
**Points:** 100

---

## Overview

Use simulation to understand the multiple testing problem and evaluate correction methods.

---

## Part A: The Problem (30 points)

### Task 1: Simulate Null Data

Generate expression data where NO genes are truly differentially expressed.

```r
set.seed(42)
n_genes <- 10000
n_samples <- 6  # 3 per group

# Generate random counts (no true differences)
counts <- matrix(
  rnbinom(n_genes * n_samples, mu = 100, size = 10),
  nrow = n_genes
)

# Assign groups
group <- factor(c(rep("A", 3), rep("B", 3)))
```

### Task 2: Run t-tests

For each gene, run a t-test comparing groups A and B.

**Questions:**
1. How many genes have p < 0.05?
2. How many genes have p < 0.01?
3. What proportion is this out of 10,000?
4. How does this compare to α (the expected false positive rate)?

---

## Part B: Multiple Testing Corrections (40 points)

### Task 3: Apply Corrections

Apply the following corrections to your p-values:

1. **Bonferroni correction**
```r
p_bonf <- p.adjust(pvalues, method = "bonferroni")
```

2. **Benjamini-Hochberg (FDR)**
```r
p_fdr <- p.adjust(pvalues, method = "BH")
```

**Questions:**
1. How many genes are significant (adjusted p < 0.05) after Bonferroni?
2. How many genes are significant after BH?
3. Which method is more conservative? Why?
4. What does FDR = 0.05 actually control?

### Task 4: Visualize Corrections

Create a plot showing:
- Raw p-values (histogram)
- Bonferroni-adjusted p-values
- BH-adjusted p-values

---

## Part C: True Positives (30 points)

### Task 5: Add True DEGs

Modify your simulation to include 500 truly differentially expressed genes.

```r
# Add fold change to first 500 genes in group B
# Your code here
```

### Task 6: Calculate Rates

After applying each correction method:

| Metric | No Correction | Bonferroni | BH |
|--------|---------------|------------|-----|
| True Positives | | | |
| False Positives | | | |
| False Negatives | | | |
| Sensitivity | | | |
| Specificity | | | |
| FDR (observed) | | | |

**Questions:**
1. Which method has higher sensitivity?
2. Which method has lower FDR?
3. What is the trade-off between these methods?
4. Why is BH preferred in practice?

---

## Deliverables

1. **R script or Rmd** with all code
2. **Figures:**
   - P-value histogram (null data)
   - Comparison of correction methods
3. **Written answers** to all questions

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Correct simulation setup | 15 |
| Proper t-test implementation | 15 |
| Correct application of corrections | 20 |
| Visualization quality | 15 |
| Interpretation of results | 20 |
| Code quality | 15 |

---

## Expected Insights

By the end, you should understand:
- Testing 10,000 genes at α=0.05 → ~500 false positives expected
- Bonferroni is very conservative (few false positives, many missed true DEGs)
- BH controls FDR while maintaining better power

