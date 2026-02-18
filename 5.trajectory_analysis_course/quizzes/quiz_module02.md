# Quiz: Module 2 - Pseudotime Concepts

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What is pseudotime?

- A) Actual experimental time
- B) Computational ordering of cells along a biological process
- C) Processing time for analysis
- D) Cell cycle duration

**Your answer:** _____

---

### Question 2
Pseudotime analysis is appropriate when:

- A) Cells exist in discrete, unrelated states
- B) Cells progress through a continuous biological process
- C) All cells are identical
- D) You have time-series data only

**Your answer:** _____

---

### Question 3
The "root" cell in trajectory analysis represents:

- A) The largest cell
- B) The starting point of the biological process
- C) A randomly selected cell
- D) The cell with highest expression

**Your answer:** _____

---

### Question 4
Branching trajectories indicate:

- A) Technical errors
- B) Multiple possible cell fates / differentiation paths
- C) Poor data quality
- D) Incorrect clustering

**Your answer:** _____

---

### Question 5
A key assumption of pseudotime analysis is:

- A) All cells were collected at the same time
- B) The snapshot captures cells at different stages of a process
- C) Cells don't change over time
- D) Gene expression is constant

**Your answer:** _____

---

### Question 6
Why can't we use actual time in single-cell experiments?

- A) Cells are destroyed during sequencing (snapshot data)
- B) Time is not important
- C) Sequencing is instantaneous
- D) All cells are synchronized

**Your answer:** _____

---

### Question 7
Pseudotime values typically range from:

- A) -100 to 100
- B) 0 to 1 (or 0 to max, normalized)
- C) Only integers
- D) Negative values only

**Your answer:** _____

---

### Question 8
What biological processes are suitable for trajectory analysis?

- A) Static tissue samples
- B) Differentiation, development, cell cycle, activation
- C) Only embryonic development
- D) Only cancer

**Your answer:** _____

---

### Question 9
Selecting the wrong root cell results in:

- A) No trajectory
- B) Reversed or incorrect pseudotime ordering
- C) Better results
- D) Faster computation

**Your answer:** _____

---

### Question 10
Validating pseudotime ordering can be done by:

- A) Checking known early/late marker genes
- B) Measuring cell size
- C) Counting clusters
- D) Ignoring biology

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
In one sentence, explain why single-cell data is "snapshot" data and what that implies for trajectory analysis.

**Your answer:**

---

### Question 12
How would you choose the root cell for a differentiation trajectory?

**Your answer:**

---

### Question 13
What does a "linear" vs "branching" trajectory tell you about the biological process?

**Your answer:**

---

### Question 14
Name two methods for inferring pseudotime (e.g., DPT, Monocle, Slingshot).

**Your answer:**

---

### Question 15
Why is root selection critical, and what happens if it is wrong?

**Your answer:**

---

## Answer Key

### Part A
1. B  
2. B  
3. B  
4. B  
5. B  
6. A  
7. B  
8. B  
9. B  
10. A  

### Part B (sample answers)
11. Cells are lysed during sequencing; we get one timepoint per cell. We infer ordering from expression similarity.  
12. Use prior knowledge (e.g. stem cell markers), or pick the cell type known to be earliest (e.g. progenitor).  
13. Linear: one path. Branching: multiple fates; cells can commit to different endpoints.  
14. Diffusion pseudotime (DPT), Monocle 3, Slingshot, destiny, Palantir.  
15. Root defines the start; wrong root flips or distorts the ordering; validation with known early/late markers is essential.
