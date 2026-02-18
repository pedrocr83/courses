# Quiz: Module 5 - Hypothesis Testing Fundamentals

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What does the null hypothesis (H₀) typically state in DE analysis?

- A) The gene is differentially expressed
- B) There is no difference in expression between groups
- C) The p-value is less than 0.05
- D) The effect size is large

**Your answer:** _____

---

### Question 2
A p-value of 0.03 means:

- A) There is a 3% chance the null hypothesis is true
- B) 3% of genes are differentially expressed
- C) If H₀ is true, there's a 3% chance of seeing data this extreme
- D) The effect size is 0.03

**Your answer:** _____

---

### Question 3
What is a Type I error?

- A) Failing to detect a true difference
- B) Rejecting H₀ when it is actually true (false positive)
- C) A computational error
- D) Using the wrong statistical test

**Your answer:** _____

---

### Question 4
What is a Type II error?

- A) Failing to reject H₀ when it is false (false negative)
- B) Rejecting H₀ when it is true
- C) Using too many samples
- D) Incorrect normalization

**Your answer:** _____

---

### Question 5
Statistical power is:

- A) The probability of a Type I error
- B) The probability of detecting a true effect
- C) The p-value threshold
- D) The fold change magnitude

**Your answer:** _____

---

### Question 6
Which test is commonly used for comparing two groups of expression values?

- A) Chi-square test
- B) t-test or Wilcoxon rank-sum test
- C) ANOVA only
- D) Correlation test

**Your answer:** _____

---

### Question 7
A parametric test (like t-test) assumes:

- A) Data is normally distributed
- B) Data has no distribution
- C) Only integer values
- D) Equal sample sizes always

**Your answer:** _____

---

### Question 8
When is a Wilcoxon test preferred over a t-test?

- A) When data is normally distributed
- B) When data violates normality assumptions
- C) When sample size is very large
- D) When comparing more than 2 groups

**Your answer:** _____

---

### Question 9
What does a very small p-value (e.g., 1e-50) indicate?

- A) A biologically important gene
- B) Strong statistical evidence against H₀
- C) The largest effect size
- D) No need for multiple testing correction

**Your answer:** _____

---

### Question 10
Why is the p-value alone insufficient for assessing biological relevance?

- A) It doesn't account for effect size
- B) It's always too large
- C) It requires normal distribution
- D) It only works for bulk RNA-seq

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
Explain the difference between Type I and Type II errors in one sentence each.

**Your answer:**

---

### Question 12
What is the correct interpretation of a p-value? (What does p=0.01 mean?)

**Your answer:**

---

### Question 13
Why might a gene with p=0.001 and log2FC=0.1 be less biologically interesting than one with p=0.02 and log2FC=2?

**Your answer:**

---

### Question 14
What assumption does the Wilcoxon rank-sum test relax compared to the t-test?

**Your answer:**

---

### Question 15
How does sample size affect statistical power?

**Your answer:**

---

## Answer Key

### Part A
1. B  
2. C  
3. B  
4. A  
5. B  
6. B  
7. A  
8. B  
9. B  
10. A  

### Part B (sample answers)
11. Type I: false positive – reject H₀ when it's true. Type II: false negative – fail to reject H₀ when it's false.  
12. If the null were true, we'd see data this extreme (or more) 1% of the time. It is NOT the probability the null is true.  
13. Effect size matters. log2FC=0.1 is tiny biologically; p=0.02 with large fold change may be more relevant.  
14. Wilcoxon doesn't assume normality; it uses ranks, so it's robust to non-normal and skewed distributions.  
15. Larger sample size → more power to detect true effects; smaller samples have higher Type II error rate.
