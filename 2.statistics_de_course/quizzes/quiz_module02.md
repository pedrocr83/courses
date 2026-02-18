# Quiz: Module 2 - Probability Distributions for Count Data

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What type of data does RNA-seq produce?

- A) Continuous measurements
- B) Non-negative integer counts
- C) Percentage values
- D) Binary presence/absence

**Your answer:** _____

---

### Question 2
What does "overdispersion" mean in the context of RNA-seq?

- A) Too many genes are expressed
- B) Variance exceeds the mean
- C) Data is normally distributed
- D) Samples are too similar

**Your answer:** _____

---

### Question 3
Why is the Poisson distribution inadequate for RNA-seq data?

- A) It only works for continuous data
- B) It assumes mean equals variance, which is rarely true
- C) It cannot handle negative values
- D) It is too computationally expensive

**Your answer:** _____

---

### Question 4
Which distribution is the backbone of DESeq2 and edgeR?

- A) Normal distribution
- B) Poisson distribution
- C) Negative binomial distribution
- D) Uniform distribution

**Your answer:** _____

---

### Question 5
What additional parameter does the negative binomial have compared to Poisson?

- A) Mean parameter
- B) Dispersion parameter
- C) Location parameter
- D) Scale parameter

**Your answer:** _____

---

### Question 6
Zero-inflated models are particularly relevant for:

- A) Bulk RNA-seq
- B) Single-cell RNA-seq
- C) Microarray data
- D) Proteomics

**Your answer:** _____

---

### Question 7
If a gene has mean expression of 100 and variance of 500, this suggests:

- A) Underdispersion
- B) Perfect Poisson behavior
- C) Overdispersion
- D) Normal distribution

**Your answer:** _____

---

### Question 8
What is the relationship between mean and variance in a Poisson distribution?

- A) Variance = Mean²
- B) Variance = Mean
- C) Variance = Mean / 2
- D) No relationship

**Your answer:** _____

---

### Question 9
In the negative binomial, higher dispersion means:

- A) Lower variance
- B) Higher variance relative to mean
- C) More Poisson-like behavior
- D) Smaller sample size needed

**Your answer:** _____

---

### Question 10
Why do RNA-seq count data typically show overdispersion?

- A) Sequencing errors only
- B) Biological variation between samples
- C) Computational artifacts
- D) Poor experimental design

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
In one sentence, explain why biological variation between replicates leads to overdispersion in RNA-seq counts.

**Your answer:**

---

### Question 12
What does the dispersion parameter in the negative binomial control, and how does it relate to the mean-variance relationship?

**Your answer:**

---

### Question 13
Why are zero-inflated models often used for single-cell but less for bulk RNA-seq?

**Your answer:**

---

### Question 14
If variance = mean for a gene, which distribution would fit best?

**Your answer:**

---

### Question 15
Name one reason DESeq2 and edgeR use the negative binomial rather than Poisson for RNA-seq.

**Your answer:**

---

## Answer Key

### Part A
1. B  
2. B  
3. B  
4. C  
5. B  
6. B  
7. C  
8. B  
9. B  
10. B  

### Part B (sample answers)
11. Biological variation means replicates differ in true expression; counts vary more than expected under Poisson (which assumes a single true rate).  
12. Dispersion controls how much variance exceeds the mean; higher dispersion = more overdispersion, variance further above mean.  
13. Single-cell has many zeros (dropouts, low capture); zero-inflated models explicitly model excess zeros. Bulk has fewer zeros.  
14. Poisson.  
15. RNA-seq data are overdispersed (variance > mean); negative binomial models this; Poisson would underestimate uncertainty.
