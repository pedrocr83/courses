# Quiz: Module 4 - CellPhoneDB

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What is the primary statistical method used by CellPhoneDB?

- A) t-test
- B) Permutation testing
- C) Chi-square test
- D) Linear regression

**Your answer:** _____

---

### Question 2
CellPhoneDB requires input data to be:

- A) Raw counts only
- B) Normalized expression + cell type annotations
- C) Protein abundance data
- D) Spatial coordinates

**Your answer:** _____

---

### Question 3
What does a low p-value in CellPhoneDB output indicate?

- A) The interaction is not significant
- B) The L-R pair is expressed higher in specific cell type pairs than expected by chance
- C) The genes are not expressed
- D) The database entry is incorrect

**Your answer:** _____

---

### Question 4
CellPhoneDB's database primarily contains:

- A) Only secreted ligands
- B) Curated ligand-receptor and receptor-receptor interactions
- C) Transcription factor networks
- D) Metabolic pathways

**Your answer:** _____

---

### Question 5
The "means" value in CellPhoneDB output represents:

- A) Average p-value
- B) Average expression of L-R pair across cell type combination
- C) Number of cells expressing
- D) Fold change

**Your answer:** _____

---

### Question 6
How does CellPhoneDB handle multi-subunit receptors?

- A) Ignores them
- B) Requires all subunits to be expressed
- C) Uses only one subunit
- D) Averages all subunits

**Your answer:** _____

---

### Question 7
Which language is CellPhoneDB primarily written in?

- A) R
- B) Python
- C) Julia
- D) Perl

**Your answer:** _____

---

### Question 8
What is a limitation of CellPhoneDB's permutation approach?

- A) It's too slow
- B) Shuffling cell labels may not preserve cluster structure
- C) Only works for mouse data
- D) Requires spatial data

**Your answer:** _____

---

### Question 9
CellPhoneDB's dot plot typically shows:

- A) Gene expression only
- B) Interaction significance and expression magnitude
- C) Cell type frequencies
- D) Spatial distances

**Your answer:** _____

---

### Question 10
When comparing conditions in CellPhoneDB, you should:

- A) Run once on combined data
- B) Run separately on each condition and compare
- C) Only analyze one condition
- D) Ignore condition information

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
How does CellPhoneDB's permutation test work? What null hypothesis does it test?

**Your answer:**

---

### Question 12
Why does CellPhoneDB require cell type annotations as input?

**Your answer:**

---

### Question 13
What is the difference between the "means" and "pvalues" columns in CellPhoneDB output?

**Your answer:**

---

### Question 14
What does it mean if an L-R pair has low p-value but low "means"?

**Your answer:**

---

### Question 15
When would you filter CellPhoneDB results by both p-value and expression (means)?

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
7. B  
8. B  
9. B  
10. B  

### Part B (sample answers)
11. Permutes cell type labels; tests whether observed L-R expression pattern is more extreme than random; null: no association between cell types and L-R expression.  
12. To define sender and receiver populations; CCC is between cell types.  
13. means: average expression level. pvalues: statistical significance of the interaction.  
14. Significant by chance structure but low actual expression; may be biologically irrelevant.  
15. To get high-confidence, high-expression interactions; filter padj < 0.05 and means > threshold.
