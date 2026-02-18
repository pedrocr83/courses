# Quiz: Module 5 - Graph-Based Clustering

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What does k-NN stand for in the context of clustering?

- A) k-Normalized Networks
- B) k-Nearest Neighbors
- C) Kernel Neural Networks
- D) k-Node Networks

**Your answer:** _____

---

### Question 2
The Leiden algorithm improves upon Louvain by:

- A) Being faster
- B) Guaranteeing connected communities
- C) Using fewer parameters
- D) Only working on small datasets

**Your answer:** _____

---

### Question 3
What does the "resolution" parameter control in Leiden clustering?

- A) Image quality of plots
- B) Number/size of clusters (higher = more clusters)
- C) Running time
- D) Memory usage

**Your answer:** _____

---

### Question 4
An SNN (Shared Nearest Neighbor) graph differs from k-NN by:

- A) Using fewer neighbors
- B) Weighting edges by shared neighbor overlap
- C) Being undirected
- D) Requiring labels

**Your answer:** _____

---

### Question 5
If you set resolution too high, you will likely:

- A) Merge distinct cell types
- B) Over-cluster, splitting real populations
- C) Get exactly the right number of clusters
- D) Crash the algorithm

**Your answer:** _____

---

### Question 6
If you set resolution too low, you will likely:

- A) Over-cluster rare populations
- B) Under-cluster, merging distinct cell types
- C) Get too many clusters
- D) Improve biological accuracy

**Your answer:** _____

---

### Question 7
The clustering input in Scanpy/Seurat is typically:

- A) Raw counts
- B) PCA embedding
- C) Gene names
- D) Cell annotations

**Your answer:** _____

---

### Question 8
What is modularity in community detection?

- A) Code organization
- B) Measure of within-cluster vs between-cluster connections
- C) Number of modules loaded
- D) GPU utilization

**Your answer:** _____

---

### Question 9
A good strategy for choosing resolution is:

- A) Always use 1.0
- B) Test multiple values and evaluate biologically
- C) Use the maximum possible value
- D) Let the algorithm decide automatically

**Your answer:** _____

---

### Question 10
The "clustree" package helps you:

- A) Cluster trees (plants)
- B) Visualize how clusters change across resolutions
- C) Speed up clustering
- D) Annotate cell types

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
In one sentence, what is the purpose of building a k-NN or SNN graph before clustering?

**Your answer:**

---

### Question 12
Why does Leiden guarantee well-connected communities whereas Louvain may not?

**Your answer:**

---

### Question 13
What is the trade-off between using too many vs too few neighbors when building the graph?

**Your answer:**

---

### Question 14
How would you decide if a cluster should be sub-clustered?

**Your answer:**

---

### Question 15
What metric could you use to evaluate clustering quality objectively (e.g., silhouette)?

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
11. To capture local similarity structure so clustering can find groups of cells that are similar to each other.  
12. Leiden uses a refinement step that ensures each community is locally connected; Louvain can produce disconnected communities.  
13. Few neighbors: sparse graph, may miss structure. Many neighbors: dense graph, may merge distinct populations.  
14. If it expresses markers of multiple cell types, appears heterogeneous on UMAP, or has high internal diversity.  
15. Silhouette score, modularity, or biological validation (known markers, expected cell types).
