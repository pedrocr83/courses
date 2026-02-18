# Quiz: Module 7 - LIANA Consensus Framework

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What problem does LIANA address in CCC analysis?

- A) Lack of databases
- B) Inconsistent results across different CCC methods
- C) Slow computation
- D) Missing cell types

**Your answer:** _____

---

### Question 2
LIANA provides a unified interface to run:

- A) Only CellPhoneDB
- B) Multiple CCC methods with consistent input/output
- C) Only R-based tools
- D) Spatial analysis only

**Your answer:** _____

---

### Question 3
LIANA's aggregate score combines:

- A) Only p-values
- B) Rankings from multiple methods
- C) Cell counts
- D) Gene lengths

**Your answer:** _____

---

### Question 4
Which database resource does LIANA access by default?

- A) Only CellChatDB
- B) OmniPath and multiple curated resources
- C) KEGG only
- D) Custom user databases only

**Your answer:** _____

---

### Question 5
Why is consensus scoring useful?

- A) It's faster
- B) Interactions identified by multiple methods are more reliable
- C) It uses less memory
- D) It requires fewer cells

**Your answer:** _____

---

### Question 6
LIANA is available in which programming languages?

- A) R only
- B) Python only
- C) Both R and Python
- D) Julia only

**Your answer:** _____

---

### Question 7
When an interaction is ranked highly by LIANA but not by individual tools, this suggests:

- A) The interaction is definitely real
- B) It may be an artifact of the consensus method
- C) The databases are incomplete
- D) The data is poor quality

**Your answer:** _____

---

### Question 8
LIANA's "select_resource" function allows you to:

- A) Download new data
- B) Choose which L-R database to use
- C) Select cell types
- D) Filter genes

**Your answer:** _____

---

### Question 9
Compared to running individual tools, LIANA:

- A) Takes much longer
- B) Provides standardized output for easier comparison
- C) Uses completely different algorithms
- D) Only works on mouse data

**Your answer:** _____

---

### Question 10
LIANA's robust rank aggregate is based on:

- A) Simple averaging
- B) Rank aggregation statistical methods
- C) Neural networks
- D) Random sampling

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
In one sentence, explain why different CCC tools often give different results.

**Your answer:**

---

### Question 12
What does "consensus" mean in the context of LIANA, and how does it improve confidence?

**Your answer:**

---

### Question 13
When might you prioritize a method-specific rank over the LIANA aggregate score?

**Your answer:**

---

### Question 14
Name two CCC methods that LIANA can run.

**Your answer:**

---

### Question 15
What is the benefit of LIANA's standardized output format across methods?

**Your answer:**

---

## Answer Key

### Part A
1. B  
2. B  
3. B  
4. B  
5. B  
6. C  
7. B  
8. B  
9. B  
10. B  

### Part B (sample answers)
11. Different scoring (permutation, mass action, etc.), databases, and null models produce different rankings.  
12. Consensus = agreement across methods; interactions ranked highly by multiple methods are more likely to be real.  
13. When one method is more appropriate for your data type; when exploring method-specific biology.  
14. CellPhoneDB, CellChat, connectome, NATMI, SingleCellSignalR, etc.  
15. Direct comparison of method outputs; same columns/format; easier to merge and analyze.
