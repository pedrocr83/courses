# Assignment 4: Build and Analyze a CCC Network

**Module:** 10 - CCC as a Graph Problem  
**Due:** End of Week 4  
**Points:** 100

---

## Overview

Construct a cell-cell communication network from your analysis and apply graph-theoretic metrics to identify key communication hubs.

---

## Dataset

Use CCC results from previous assignments (CellChat or LIANA output).

---

## Part A: Network Construction (30 points)

### Task 1: Define Network Structure (15 points)

Choose one of two approaches:

**Option 1: Cell Type Network**
- Nodes = Cell types
- Edges = Aggregate interactions between cell types
- Edge weight = Number of interactions or sum of scores

**Option 2: L-R Interaction Network**
- Nodes = Ligand-receptor pairs
- Edges = Shared cell type involvement
- More complex but more detailed

Document your choice and rationale.

### Task 2: Build the Network (15 points)

```r
library(igraph)

# From CellChat
net <- cellchat@net$weight  # Cell type × Cell type matrix

# Create igraph object
g <- graph_from_adjacency_matrix(
  net, 
  mode = "directed", 
  weighted = TRUE
)
```

Or from a generic interaction table:

```r
# From edge list
edges <- data.frame(
  from = interaction_table$sender,
  to = interaction_table$receiver,
  weight = interaction_table$score
)

g <- graph_from_data_frame(edges, directed = TRUE)
```

---

## Part B: Network Metrics (35 points)

### Task 3: Calculate Centrality Measures (20 points)

For each cell type (node):

1. **Out-degree:** Number of outgoing interactions (sender activity)
2. **In-degree:** Number of incoming interactions (receiver activity)
3. **Betweenness:** How often a node bridges others
4. **PageRank:** Importance based on connection quality

```r
# Calculate metrics
V(g)$out_degree <- degree(g, mode = "out")
V(g)$in_degree <- degree(g, mode = "in")
V(g)$betweenness <- betweenness(g)
V(g)$pagerank <- page_rank(g)$vector
```

### Task 4: Identify Key Nodes (15 points)

1. Which cell type is the top sender?
2. Which is the top receiver?
3. Which has highest betweenness (bridge/hub)?
4. Do these align with CellChat's signaling role analysis?

---

## Part C: Visualization (20 points)

### Task 5: Create Network Visualizations (20 points)

1. **Force-directed layout** with:
   - Node size = degree or PageRank
   - Edge width = interaction strength
   - Node color = cell type

2. **Circular layout** (similar to CellChat)

3. **Heatmap** of the adjacency matrix

```r
# Basic plot
plot(g, 
     vertex.size = V(g)$out_degree * 2,
     edge.width = E(g)$weight / max(E(g)$weight) * 5,
     layout = layout_with_fr(g))
```

---

## Part D: Interpretation (15 points)

### Task 6: Biological Interpretation (15 points)

Write 200-300 words addressing:

1. What does the network structure tell you about this tissue?
2. Are the hub cell types biologically expected?
3. What would happen if you removed the top hub?
4. How might this network differ in disease?

---

## Deliverables

1. **R or Python script** with network construction and analysis
2. **Figures:**
   - Network visualization (at least 2 layouts)
   - Adjacency heatmap
   - Centrality barplots
3. **Tables:**
   - Node metrics for all cell types
   - Top 5 nodes by each metric
4. **Written interpretation**

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Network construction correct | 15 |
| Network structure justified | 15 |
| Centrality metrics calculated | 20 |
| Key nodes identified | 15 |
| Visualizations quality | 20 |
| Biological interpretation | 15 |

---

## Advanced Extension (Optional)

- Detect communities in the network
- Compare networks between conditions
- Simulate signaling propagation

