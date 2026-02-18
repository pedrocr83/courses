# Quiz: Module 5 - CellChat

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
CellChat uses what type of model to score interactions?

- A) Permutation-based
- B) Probability-based using law of mass action
- C) Correlation-based
- D) Simple co-expression

**Your answer:** _____

---

### Question 2
What is "communication probability" in CellChat?

- A) p-value from a statistical test
- B) A score combining expression and interaction strength
- C) Percentage of cells expressing
- D) Fold change between conditions

**Your answer:** _____

---

### Question 3
CellChat aggregates interactions into:

- A) Individual gene pairs only
- B) Signaling pathways
- C) Chromosomal regions
- D) Metabolic modules

**Your answer:** _____

---

### Question 4
Which cell roles does CellChat identify in signaling networks?

- A) Senders, receivers, mediators, influencers
- B) Only senders and receivers
- C) Only hubs
- D) Only leaf nodes

**Your answer:** _____

---

### Question 5
CellChat's circle plot shows:

- A) Gene expression levels
- B) Direction and strength of communication between cell types
- C) Cell clustering
- D) Spatial locations

**Your answer:** _____

---

### Question 6
CellChatDB includes which types of interactions?

- A) Secreted signaling only
- B) Secreted, ECM-receptor, and cell-cell contact
- C) Only receptor-receptor
- D) Only intracellular

**Your answer:** _____

---

### Question 7
To compare two conditions in CellChat, you should:

- A) Combine datasets and run once
- B) Create separate CellChat objects and merge
- C) Only analyze one condition
- D) Use a different tool

**Your answer:** _____

---

### Question 8
What visualization shows hierarchical communication patterns in CellChat?

- A) Volcano plot
- B) River plot / Sankey diagram
- C) PCA plot
- D) MA plot

**Your answer:** _____

---

### Question 9
CellChat requires which R package for visualization?

- A) ggplot2 only
- B) ComplexHeatmap and circlize
- C) Shiny
- D) plotly

**Your answer:** _____

---

### Question 10
The "netAnalysis_signalingRole" function identifies:

- A) Differentially expressed genes
- B) Dominant senders and receivers based on network centrality
- C) Cell type markers
- D) Quality control metrics

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
How does CellChat's law of mass action differ from CellPhoneDB's permutation approach?

**Your answer:**

---

### Question 12
What advantage does pathway-level aggregation provide over individual L-R pairs?

**Your answer:**

---

### Question 13
What do "influencer" and "mediator" cells mean in CellChat's network analysis?

**Your answer:**

---

### Question 14
Why might CellChat and CellPhoneDB give different results for the same dataset?

**Your answer:**

---

### Question 15
When would you use CellChat's merge comparison to analyze two conditions?

**Your answer:**

---

## Answer Key

### Part A
1. B  
2. B  
3. B  
4. A  
5. B  
6. B  
7. B  
8. B  
9. B  
10. B  

### Part B (sample answers)
11. CellChat uses expression levels and interaction probabilities (mass action); CellPhoneDB uses permutation to test significance. Different statistical frameworks.  
12. Multiple L-R pairs can act through same pathway; pathway-level is more interpretable and robust.  
13. Influencer: cells that affect many others. Mediator: bridge between sender and receiver populations.  
14. Different scoring (mass action vs permutation), different databases (CellChatDB vs CellPhoneDB), different pathway aggregation.  
15. When you have control vs treatment (or condition A vs B) and want to compare which pathways/L-R are stronger in each.
