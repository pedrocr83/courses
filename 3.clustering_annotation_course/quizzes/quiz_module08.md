# Quiz: Module 8 - Manual Annotation

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What is the primary basis for manual cell type annotation?

- A) Cluster numbers
- B) Known marker gene expression
- C) Cell size
- D) Sequencing depth

**Your answer:** _____

---

### Question 2
CD3D, CD3E, and CD3G are canonical markers for:

- A) B cells
- B) T cells
- C) Monocytes
- D) NK cells

**Your answer:** _____

---

### Question 3
CD14 and LYZ are markers for:

- A) T cells
- B) B cells
- C) Monocytes/Macrophages
- D) Platelets

**Your answer:** _____

---

### Question 4
MS4A1 (CD20) is a marker for:

- A) T cells
- B) B cells
- C) Monocytes
- D) Dendritic cells

**Your answer:** _____

---

### Question 5
When a cluster expresses markers from multiple cell types, this suggests:

- A) A novel cell type
- B) Possible doublets or need for sub-clustering
- C) Perfect annotation
- D) Technical error - restart analysis

**Your answer:** _____

---

### Question 6
PanglaoDB and CellMarker are:

- A) Clustering algorithms
- B) Databases of cell type marker genes
- C) Normalization methods
- D) Programming languages

**Your answer:** _____

---

### Question 7
What is a "dot plot" useful for in annotation?

- A) Showing cluster sizes
- B) Visualizing marker expression across clusters
- C) PCA results
- D) Trajectory analysis

**Your answer:** _____

---

### Question 8
If a known T cell marker shows expression in your "B cell" cluster, you should:

- A) Ignore it
- B) Investigate - possible mis-annotation, doublets, or contamination
- C) Delete the cluster
- D) Change the marker database

**Your answer:** _____

---

### Question 9
Sub-clustering is useful when:

- A) You have too few cells
- B) A cluster appears heterogeneous with multiple subtypes
- C) The algorithm failed
- D) You want fewer clusters

**Your answer:** _____

---

### Question 10
Good annotation practice includes:

- A) Only using automated tools
- B) Documenting evidence for each annotation decision
- C) Assigning names randomly
- D) Ignoring ambiguous clusters

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
Name two visualization types used to assess marker expression for annotation.

**Your answer:**

---

### Question 12
Why should you use multiple markers rather than a single marker to assign a cell type?

**Your answer:**

---

### Question 13
What does "confidence" mean in the context of cell type annotation, and how might you assign it?

**Your answer:**

---

### Question 14
When might you label a cluster as "unknown" or "unassigned" rather than forcing an annotation?

**Your answer:**

---

### Question 15
How do automated annotation tools (e.g., SingleR, CellTypist) complement manual annotation?

**Your answer:**

---

## Answer Key

### Part A
1. B  
2. B  
3. C  
4. B  
5. B  
6. B  
7. B  
8. B  
9. B  
10. B  

### Part B (sample answers)
11. Dot plot, violin plot, feature plot, heatmap.  
12. Single markers can be shared or promiscuous; multiple markers increase confidence and specificity.  
13. Confidence reflects how strongly the evidence supports the annotation; high = strong marker support; low = ambiguous or conflicting.  
14. When markers are unclear, cluster is heterogeneous, or no good reference exists; better to be honest than wrong.  
15. They provide a starting point or validation; manual review catches errors and adds biological nuance.
