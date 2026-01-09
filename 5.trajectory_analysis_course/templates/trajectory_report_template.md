# Trajectory Analysis Report

## Dataset Overview

| Property | Value |
|----------|-------|
| Dataset | |
| Cells | |
| Cell types | |
| Biological process | |
| Expected trajectory | |

---

## Biological Context

**Process being studied:** 

**Expected stages:**
1. 
2. 
3. 

**Key genes expected to change:**
- Early: 
- Late: 

---

## Trajectory Inference

### Method 1: [Diffusion PT / Monocle / Slingshot]

| Parameter | Value |
|-----------|-------|
| Root cell type | |
| | |

**Results:**
- Pseudotime range: 
- Terminal states identified: 

### Method 2: [Name]

| Parameter | Value |
|-----------|-------|
| | |

**Results:**

### Method Comparison

| Aspect | Method 1 | Method 2 |
|--------|----------|----------|
| Correlation | | |
| Root agreement | | |
| Terminal agreement | | |

**Selected method:** 
**Rationale:** 

---

## RNA Velocity Analysis

### Velocity Computation

| Parameter | Value |
|-----------|-------|
| Model | Stochastic/Dynamical |
| Genes with velocity | |
| Mean R² | |

### Velocity Interpretation

| Cell Type | Velocity Magnitude | Direction |
|-----------|-------------------|-----------|
| | High/Low | Toward X |
| | | |

**Agreement with pseudotime:** [Yes/Partial/No]

---

## Trajectory Differential Expression

### Summary

| Category | Gene Count |
|----------|------------|
| Total significant | |
| Increasing | |
| Decreasing | |
| Transient | |

### Top Trajectory Genes

#### Increasing (late markers)
1. GENE1 - pattern, function
2. 
3. 

#### Decreasing (early markers)
1. 
2. 
3. 

#### Transient (intermediate)
1. 
2. 

---

## Gene Dynamics Visualization

[Heatmap along pseudotime]

[Top gene expression curves]

---

## CellRank / Fate Analysis

### Terminal States

| State | Cells | Fate Probability |
|-------|-------|------------------|
| | | |
| | | |

### Driver Genes

| Terminal State | Top Drivers |
|----------------|-------------|
| | |
| | |

---

## Pathway Enrichment

### Genes increasing along trajectory

| Pathway | p-value | Genes |
|---------|---------|-------|
| | | |
| | | |

### Genes decreasing

| Pathway | p-value | Genes |
|---------|---------|-------|
| | | |

---

## Biological Interpretation

**Key findings:**
1. 
2. 
3. 

**Agreement with literature:**

**Novel observations:**

---

## Limitations

| Limitation | Impact |
|------------|--------|
| | |
| | |

---

## Files Generated

- `adata_trajectory.h5ad`
- `pseudotime.csv`
- `trajectory_de_genes.csv`
- `velocity_genes.csv`
- Figures in `/figures/`

