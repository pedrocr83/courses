# Module 11A: Pseudobulk Theory & When to Use

**Understanding the statistical foundation of pseudobulk analysis for scRNA-seq**

---

## Overview

Pseudobulk analysis is **the recommended approach** for differential expression (DE) analysis in single-cell RNA-seq when comparing across samples/conditions. This module explains why, when, and how to properly implement pseudobulk methods.

**Duration:** 1.5-2 hours  
**Prerequisites:** Modules 1-8 (bulk RNA-seq DE), basic understanding of scRNA-seq

---

## Learning Objectives

By the end of this module, you will be able to:
- [ ] Explain why cell-level DE is statistically problematic
- [ ] Understand the concept of pseudoreplication
- [ ] Correctly aggregate cells to create pseudobulk samples
- [ ] Choose appropriate aggregation strategies
- [ ] Apply bulk DE methods (DESeq2) to pseudobulk data
- [ ] Interpret pseudobulk DE results correctly

---

## The Problem with Cell-Level DE

### What is Pseudoreplication?

**Pseudoreplication** occurs when technical replicates (cells from the same sample) are treated as independent biological replicates.

#### Example Scenario

**Experiment:**
- 3 Control mice
- 3 Treated mice
- 10,000 cells per mouse (60,000 cells total)

**WRONG Analysis (Cell-Level):**
```python
# Treating each cell as independent
# Group 1: 30,000 control cells
# Group 2: 30,000 treated cells
# n = 30,000 per group ❌ WRONG!
```

**Problem:** Cells from the same mouse are **not independent**. They share:
- Genetic background
- Environmental exposures
- Sample processing conditions
- Capture/sequencing batch
- Biological donor variance

**Consequence:** Massively inflated sample size → artificially tiny p-values → false discoveries

---

### Mathematical Demonstration

#### True vs Pseudoreplicated Sample Size

**Truth:**
- **n = 3** (biological replicates = mice)
- Degrees of freedom = 3 + 3 - 2 = 4

**Cell-level (pseudoreplication):**
- **n = 30,000** (cells treated as replicates)
- Degrees of freedom = 60,000 - 2 = 59,998

**Impact on p-values:**

```python
# Simulation
import numpy as np
from scipy import stats

# True biological variance (between mice)
mouse_means = np.random.normal(0, 1, size=3)  # 3 mice, SD=1
# Cells from each mouse
control_cells = np.concatenate([np.random.normal(m, 0.2, 10000) for m in mouse_means])

mouse_means_treat = np.random.normal(0.1, 1, size=3)  # Slight shift
treated_cells = np.concatenate([np.random.normal(m, 0.2, 10000) for m in mouse_means_treat])

# WRONG: t-test on cells
t_cells, p_cells = stats.ttest_ind(control_cells, treated_cells)
print(f"Cell-level p-value: {p_cells:.2e}")  # Likely < 1e-10

# CORRECT: t-test on mouse means
control_means = [np.mean(np.random.normal(m, 0.2, 10000)) for m in mouse_means]
treated_means = [np.mean(np.random.normal(m, 0.2, 10000)) for m in mouse_means_treat]
t_mouse, p_mouse = stats.ttest_ind(control_means, treated_means)
print(f"Mouse-level p-value: {p_mouse:.2f}")  # Much larger, appropriate
```

**Result:**
- Cell-level: p ≈ 10⁻¹⁰ (falsely "significant")
- Mouse-level: p ≈ 0.50 (correct, not significant)

---

### Why This Matters

**Real Example:**
- Study: Mouse model of disease
- 2 diseased mice, 2 control mice
- 5,000 T cells per mouse

**Cell-level DE:**
- Claims 5,000 differentially expressed genes (FDR < 0.05)
- But you only have **n=2 per group**!
- Impossible to detect biological variability with n=2

**Pseudobulk DE:**
- Claims 100 differentially expressed genes (FDR < 0.05)
- More realistic given limited statistical power
- Properly accounts for donor-to-donor variation

---

## The Pseudobulk Solution

### Core Concept

**Aggregate cells within each sample, then apply bulk DE methods**

```
Sample 1 (Control, Mouse 1):  10,000 cells → 1 pseudobulk profile
Sample 2 (Control, Mouse 2):  10,000 cells → 1 pseudobulk profile
Sample 3 (Control, Mouse 3):  10,000 cells → 1 pseudobulk profile
Sample 4 (Treated, Mouse 4):  10,000 cells → 1 pseudobulk profile
Sample 5 (Treated, Mouse 5):  10,000 cells → 1 pseudobulk profile
Sample 6 (Treated, Mouse 6):  10,000 cells → 1 pseudobulk profile

→ 6 pseudobulk samples for DESeq2
```

**Key insight:** The **sample** (mouse, donor, biological replicate) is the experimental unit, not the cell.

---

### When to Use Pseudobulk

✅ **USE PSEUDOBULK WHEN:**
- Comparing conditions across biological replicates (samples/donors)
- You have ≥3 samples per condition
- Research question is about population-level changes
- You want to generalize beyond your specific samples

Examples:
- Disease vs healthy (across patients)
- Treated vs untreated (across mice)
- Timepoint comparisons (across individuals)
- Genotype comparisons (across animals)

---

❌ **DON'T USE PSEUDOBULK WHEN:**
- Comparing cell types within the same sample (use cell-level methods)
- You have only 1-2 samples per condition (not enough power regardless)
- Research question is about cell-cell heterogeneity (use distribution-based methods)
- Comparing conditions that exist within single samples (e.g., tumor vs normal regions)

---

### Statistical Framework

#### Hierarchical Structure of scRNA-seq Data

```
Level 1: Condition (Treatment, Disease, etc.)
  └─ Level 2: Sample (Mouse, Patient, Donor)
      └─ Level 3: Cell (Individual cells within sample)
```

**Inference goal:** Generalize from **samples** to **condition** (population-level)

**Correct approach:**
1. Aggregate cells → sample-level summaries
2. Model sample-level variation
3. Test condition effect accounting for sample variance

**Mathematical model:**
```
Gene expression ~ Condition + Sample(Condition) + Cell(Sample)
```

Pseudobulk collapses `Cell(Sample)` by aggregation, leaving:
```
Pseudobulk expression ~ Condition
```

with sample-to-sample variation properly modeled.

---

## Pseudobulk Aggregation Strategies

### Strategy 1: Sum Counts (Recommended)

**Method:** Sum raw counts across all cells within sample × cell type

```python
import scanpy as sc
import pandas as pd

# For each sample and cell type, sum counts
def create_pseudobulk(adata, sample_col='sample', celltype_col='celltype'):
    """
    Aggregate cells to pseudobulk samples
    """
    pseudobulk_dict = {}
    
    for sample in adata.obs[sample_col].unique():
        for celltype in adata.obs[celltype_col].unique():
            # Subset to sample × cell type
            mask = (adata.obs[sample_col] == sample) & \
                   (adata.obs[celltype_col] == celltype)
            cells = adata[mask]
            
            if cells.n_obs > 0:  # If cells exist
                # Sum counts
                pseudobulk = cells.X.sum(axis=0).A1  # Convert sparse to dense
                pseudobulk_dict[f"{sample}_{celltype}"] = pseudobulk
    
    # Create pseudobulk matrix
    pseudobulk_df = pd.DataFrame(pseudobulk_dict, index=adata.var_names).T
    
    return pseudobulk_df

# Usage
pseudobulk = create_pseudobulk(adata, sample_col='donor', celltype_col='celltype')
```

**Why sum?**
- Preserves count nature (Poisson/NB distribution)
- Compatible with DESeq2/edgeR
- Doesn't artificially inflate/deflate variance
- Matches bulk RNA-seq data structure

---

### Strategy 2: Mean/Median (NOT Recommended for DE)

**Method:** Take mean or median expression per gene

❌ **Problems:**
- Loses count distribution (no longer integers)
- Variance structure changed
- Incompatible with count-based DE methods
- Information loss (especially with sparse data)

**When to use:** Only for visualization or exploratory analysis, NOT for DE testing

---

### Strategy 3: Weighted Aggregation

**Method:** Weight cells by quality metrics before summing

```python
def create_weighted_pseudobulk(adata, sample_col='sample', celltype_col='celltype',
                                 weight_col='total_counts'):
    """
    Aggregate with cell quality weighting
    """
    pseudobulk_dict = {}
    
    for sample in adata.obs[sample_col].unique():
        for celltype in adata.obs[celltype_col].unique():
            mask = (adata.obs[sample_col] == sample) & \
                   (adata.obs[celltype_col] == celltype)
            cells = adata[mask]
            
            if cells.n_obs > 0:
                # Weight by cell quality
                weights = cells.obs[weight_col].values
                weights = weights / weights.sum()  # Normalize
                
                # Weighted sum
                pseudobulk = (cells.X.T @ weights).A1
                pseudobulk_dict[f"{sample}_{celltype}"] = pseudobulk
    
    return pd.DataFrame(pseudobulk_dict, index=adata.var_names).T
```

**When to use:** If cell quality varies substantially within samples (use cautiously)

---

### Cell Type-Specific Pseudobulk

**Critical:** Aggregate separately for each cell type!

**Why?**
- Different cell types have different biology
- Mixing cell types dilutes signal
- Cell type composition differs across samples

**Example:**

```python
# WRONG: Aggregate all cells together
pseudobulk_all = adata.groupby('sample').sum()  # ❌ Mixes cell types

# CORRECT: Aggregate per cell type
for celltype in adata.obs['celltype'].unique():
    cells = adata[adata.obs['celltype'] == celltype]
    pseudobulk_ct = cells.groupby('sample').sum()  # ✅ Cell type-specific
    # Run DESeq2 on pseudobulk_ct
```

---

## Minimum Sample Size Requirements

### Power Analysis for Pseudobulk

**Rule of thumb:**
- **Minimum:** 3 samples per condition
- **Adequate:** 5-6 samples per condition
- **Good power:** 8-10 samples per condition
- **High power:** 12+ samples per condition

**Why?**
- Need to estimate biological variance
- With n=2, variance estimate unreliable
- DESeq2/edgeR shrinkage helps, but n≥3 is critical

### What if You Have Fewer Samples?

**n = 2 per group:**
- Very low power
- Use but report limitations
- Focus on genes with large effect sizes
- Validate findings in independent cohort

**n = 1 per group:**
- Cannot estimate variance
- Pseudobulk not applicable
- Consider cell-level methods (with caveats) or don't claim generalizability

---

## Common Mistakes & How to Avoid

### Mistake 1: Forgetting to Aggregate Per Cell Type

**Wrong:**
```python
# Aggregating all T cells together (mixing CD4+ and CD8+)
t_cells = adata[adata.obs['celltype'].isin(['CD4_T', 'CD8_T'])]
pseudobulk = t_cells.groupby('sample').sum()  # ❌ Mixed!
```

**Correct:**
```python
# Separate analyses for CD4+ and CD8+ T cells
cd4_pseudobulk = adata[adata.obs['celltype'] == 'CD4_T'].groupby('sample').sum()
cd8_pseudobulk = adata[adata.obs['celltype'] == 'CD8_T'].groupby('sample').sum()
```

---

### Mistake 2: Using Normalized Counts

**Wrong:**
```python
# Using log-normalized data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
pseudobulk = adata.groupby('sample').sum()  # ❌ Wrong! Log-transformed data
```

**Correct:**
```python
# Use raw counts (before normalization)
pseudobulk = adata.raw.to_adata().groupby('sample').sum()  # ✅ Raw counts
# DESeq2 will do its own normalization
```

**Why?** DESeq2/edgeR expect raw integer counts and perform their own normalization.

---

### Mistake 3: Not Handling Zero-Count Samples

**Problem:** Some samples may have zero cells for a given cell type

```python
# Example: Sample 3 has no NK cells
# Creates missing pseudobulk profile → analysis fails
```

**Solution:**
```python
# Filter: Only analyze cell types present in all samples
for celltype in celltypes:
    samples_with_cells = adata[adata.obs['celltype'] == celltype].obs['sample'].unique()
    
    if len(samples_with_cells) < n_samples:
        print(f"Skipping {celltype}: not present in all samples")
        continue
    
    # Proceed with analysis
```

**Or:** Require minimum cell count per sample
```python
min_cells = 10
# Only keep sample × cell type combinations with ≥10 cells
```

---

### Mistake 4: Ignoring Cell Count Differences

**Problem:** Wildly different cell numbers per sample can bias results

```python
# Sample 1: 50 T cells
# Sample 2: 5000 T cells
# Summed counts very different even without biological change
```

**Solutions:**

1. **Filter samples with too few cells:**
```python
min_cells = 20
samples_to_keep = cell_counts[cell_counts >= min_cells]
```

2. **Include total cell count as covariate:**
```python
# In DESeq2 design
design = ~ total_cells + condition
```

3. **Downsample to equal cells per sample:**
```python
import scanpy as sc
min_count = min([...])  # Minimum across samples
for sample in samples:
    sc.pp.subsample(adata[adata.obs['sample'] == sample], n_obs=min_count)
```

---

## Pseudobulk vs Other Methods

### Comparison Table

| Method | Use Case | Pros | Cons |
|--------|----------|------|------|
| **Pseudobulk + DESeq2** | Sample-level inference | Proper statistics, established method | Loses cell heterogeneity info |
| **Wilcoxon (Seurat)** | Within-sample comparisons | Fast, non-parametric | Pseudoreplication if comparing samples |
| **MAST** | Cell-level with covariates | Accounts for cellular detection rate | Still treats cells as independent |
| **Mixed models** | Complex designs | Proper hierarchical modeling | Computationally intensive, complex |
| **scVI DE** | Large-scale datasets | Scalable, deep learning | Black box, less interpretable |

---

### When Each Method is Appropriate

#### Use Pseudobulk When:
- ✅ Comparing conditions across samples/donors
- ✅ Have ≥3 biological replicates per group
- ✅ Want to generalize to population
- ✅ Standard experimental design (simple comparisons)

#### Use Cell-Level Methods When:
- ✅ Comparing cell types within same sample
- ✅ Identifying marker genes (no biological replicates needed)
- ✅ Exploratory analysis of cell heterogeneity
- ✅ Single-sample datasets (but can't generalize!)

#### Use Mixed Models When:
- ✅ Complex nested designs (e.g., cells within patients within clinics)
- ✅ Need to model cell-level AND sample-level effects
- ✅ Repeated measures or longitudinal designs
- ✅ Have computational resources and statistical expertise

---

## Practical Implementation

### Complete Pseudobulk Workflow

```python
import scanpy as sc
import pandas as pd
import numpy as np

# Step 1: Load QC'd single-cell data (use raw counts!)
adata = sc.read_h5ad('qc_filtered_data.h5ad')

# Ensure raw counts available
if adata.raw is None:
    raise ValueError("Need raw counts! Save before normalization.")

adata = adata.raw.to_adata()  # Use raw counts

# Step 2: Filter to cell type of interest
celltype = 'CD4_T_cells'
adata_ct = adata[adata.obs['celltype'] == celltype].copy()

print(f"Analyzing {celltype}: {adata_ct.n_obs} cells")

# Step 3: Check sample sizes
sample_sizes = adata_ct.obs.groupby(['sample', 'condition']).size()
print("\nCells per sample:")
print(sample_sizes)

# Filter samples with too few cells
min_cells_per_sample = 20
samples_to_keep = sample_sizes[sample_sizes >= min_cells_per_sample].index.get_level_values('sample')
adata_ct = adata_ct[adata_ct.obs['sample'].isin(samples_to_keep)]

print(f"\nKept {len(samples_to_keep)} samples with ≥{min_cells_per_sample} cells")

# Step 4: Create pseudobulk by summing counts per sample
pseudobulk_dict = {}
metadata_dict = {}

for sample in adata_ct.obs['sample'].unique():
    sample_cells = adata_ct[adata_ct.obs['sample'] == sample]
    
    # Sum counts
    pseudobulk_dict[sample] = sample_cells.X.sum(axis=0).A1
    
    # Store metadata
    metadata_dict[sample] = {
        'condition': sample_cells.obs['condition'].iloc[0],
        'n_cells': sample_cells.n_obs,
        'donor': sample_cells.obs['donor'].iloc[0] if 'donor' in sample_cells.obs else sample
    }

# Create count matrix (genes × samples)
pseudobulk_counts = pd.DataFrame(pseudobulk_dict, index=adata_ct.var_names)
pseudobulk_meta = pd.DataFrame(metadata_dict).T

print(f"\nPseudobulk matrix: {pseudobulk_counts.shape[0]} genes × {pseudobulk_counts.shape[1]} samples")
print(f"Conditions: {pseudobulk_meta['condition'].value_counts()}")

# Step 5: Save for DESeq2 (move to R)
pseudobulk_counts.to_csv('pseudobulk_counts.csv')
pseudobulk_meta.to_csv('pseudobulk_metadata.csv')

print("\nReady for DESeq2 analysis!")
```

---

### DESeq2 Analysis (in R)

```r
library(DESeq2)

# Load pseudobulk data
counts <- read.csv('pseudobulk_counts.csv', row.names=1)
metadata <- read.csv('pseudobulk_metadata.csv', row.names=1)

# Ensure integers
counts <- round(counts)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, 
               contrast = c("condition", "treated", "control"),
               alpha = 0.05)

# Filter by FDR and fold change
sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

cat("Significant DEGs:", nrow(sig_genes), "\n")

# Save results
write.csv(as.data.frame(res), 'pseudobulk_DE_results.csv')
```

---

## Quality Control for Pseudobulk

### QC Checklist

- [ ] **Sufficient samples:** ≥3 per condition
- [ ] **Sufficient cells per sample:** ≥20 cells (minimum)
- [ ] **Similar cell counts across samples:** Check variance
- [ ] **Samples cluster by condition, not batch:** PCA plot
- [ ] **Size factors reasonable:** Check normalization
- [ ] **Dispersion estimates sensible:** Check DESeq2 plots

### QC Plots

```r
# 1. PCA of pseudobulk samples
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

# Should see:
# - Samples cluster by condition
# - No strong batch effects

# 2. Size factors
sizeFactors(dds)
# Should be ~1, not orders of magnitude different

# 3. Dispersion plot
plotDispEsts(dds)
# Should see typical shrinkage pattern

# 4. MA plot
plotMA(res, ylim=c(-5, 5))
# Should see symmetric distribution around 0
```

---

## Advantages of Pseudobulk

✅ **Proper statistical inference**
- Correct sample size (n = number of samples, not cells)
- Appropriate p-values and FDR control
- Generalizable to population

✅ **Established methods**
- Use mature, well-validated tools (DESeq2, edgeR)
- Extensive documentation and support
- Known performance characteristics

✅ **Computational efficiency**
- Reduces data dimensionality
- Fast DE testing (6 samples vs 60,000 cells)
- Works on standard computers

✅ **Interpretability**
- Clear biological meaning (sample-level changes)
- Easy to explain to biologists
- Compatible with bulk RNA-seq literature

---

## Limitations of Pseudobulk

❌ **Loss of cell-level information**
- Cannot detect changes in expression variance
- Misses shifts in cell state distributions
- Aggregates over heterogeneity

❌ **Requires multiple samples**
- Need ≥3 biological replicates
- Cannot use with single-sample experiments
- Sample-level confounders matter

❌ **Cell type-specific**
- Must analyze each cell type separately
- Miss interactions between cell types
- Requires good cell type annotations

---

## Summary & Decision Tree

### When to Use Pseudobulk?

```
START: You have scRNA-seq data

Q1: Are you comparing across samples/donors?
  NO → Use cell-level methods (e.g., Wilcoxon for markers)
  YES → Continue

Q2: Do you have ≥3 samples per condition?
  NO → Very low power; consider alternatives or note limitations
  YES → Continue

Q3: Do you want sample-level (population) inference?
  NO → Use cell-level or distribution-based methods
  YES → **USE PSEUDOBULK** ✅

Q4: Do you need to model cell-level AND sample-level effects?
  NO → Pseudobulk sufficient
  YES → Consider mixed models (advanced)
```

---

## Key Takeaways

1. **Pseudobulk = correct statistics** for sample-level comparisons
2. **Pseudoreplication = major problem** in cell-level DE
3. **Aggregate per cell type** separately
4. **Use raw counts** for pseudobulk (not normalized)
5. **Need ≥3 samples** per condition minimum
6. **Apply established methods** (DESeq2/edgeR) after aggregation

---

## Practice Exercise

**Dataset:** Mouse PBMC, 3 control, 3 LPS-treated, multiple cell types

**Tasks:**
1. [ ] Identify appropriate cell types to analyze (≥20 cells/sample)
2. [ ] Create pseudobulk profiles for CD4 T cells
3. [ ] Check QC (cell counts, PCA)
4. [ ] Run DESeq2
5. [ ] Compare to cell-level DE (observe p-value inflation)

---

## Resources

### Key Papers
- Crowell et al. (2020) "muscat: Multi-sample multi-group scRNA-seq data analysis"
- Squair et al. (2021) "Confronting false discoveries in single-cell DE"
- Zimmerman et al. (2021) "A practical solution to pseudoreplication bias"

### Tools
- **muscat (R):** Purpose-built for pseudobulk
- **DESeq2 (R):** Standard bulk RNA-seq DE
- **edgeR (R):** Alternative bulk DE method
- **Scanpy (Python):** Data aggregation

---

**Next Steps:**
- See Module 11B for hands-on implementation
- Practice with `labs/lab11_pseudobulk.ipynb`
- Read Crowell et al. (2020) for comprehensive treatment

**Remember:** When comparing across samples, pseudobulk is not just an option—it's the statistically correct approach!
