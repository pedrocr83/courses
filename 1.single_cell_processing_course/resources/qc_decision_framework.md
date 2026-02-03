# Quality Control Decision Framework

**A systematic guide for making QC threshold decisions in single-cell RNA-seq processing**

---

## Overview

Quality control is **context-dependent**. There is no one-size-fits-all threshold. This framework helps you make informed, justifiable QC decisions based on your specific experiment.

**Key Principle:** *The goal of QC is to remove low-quality cells while preserving biological variation, including rare populations.*

---

## QC Decision Process

### Step 1: Understand Your Biological Context

Before setting any thresholds, answer these questions:

#### Sample Type Questions
- [ ] **Tissue type:** _________________
  - Brain/neurons → Higher mt% acceptable (10-20%)
  - Blood/PBMC → Lower mt% expected (5-10%)
  - Solid tumor → Expect high heterogeneity
  - Sorted cells → Narrower QC distributions expected

- [ ] **Cell state expectations:**
  - Metabolically active cells → Higher UMI counts
  - Quiescent cells → Lower gene counts acceptable
  - Differentiating cells → Variable expression patterns

- [ ] **Rare population presence:**
  - [ ] Do you expect rare cell types (< 1% of population)?
  - [ ] Are these biologically important?
  - [ ] Risk of filtering them out?

#### Technical Questions
- [ ] **Protocol used:** 10x v2 / v3 / Drop-seq / other: _____
- [ ] **Target cell number:** _____
- [ ] **Expected cells per sample:** _____
- [ ] **Sequencing depth per cell:** _____ reads
- [ ] **Fresh vs frozen tissue:** _____

---

## Step 2: Generate QC Metrics

Calculate these metrics for every cell:

### Essential Metrics
1. **n_genes_by_counts** (or nFeature_RNA): Number of genes with ≥1 count
2. **total_counts** (or nCount_RNA): Total UMI count per cell
3. **pct_counts_mt** (or percent.mt): % of counts from mitochondrial genes
4. **pct_counts_ribo**: % of counts from ribosomal genes (optional)
5. **log10_total_counts**: Log-transformed total counts (for visualization)

### Useful Derived Metrics
6. **log10_n_genes**: Log-transformed gene count
7. **mt_ratio**: Ratio of mt% to total counts (indicator of cell stress)
8. **complexity**: n_genes / total_counts (transcriptomic complexity)

---

## Step 3: Visualize Distributions

Create these essential plots:

### Plot 1: Violin Plots by Sample
```python
# Scanpy example
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, groupby='sample')
```

**Look for:**
- Bimodal distributions (may indicate cell types OR quality issues)
- Extreme outliers
- Sample-to-sample variation
- Long low-quality tails

### Plot 2: Scatter Plots (Relationships)
```python
# Total counts vs genes detected
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt')

# Mt% vs total counts
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
```

**Look for:**
- Linear relationship between total counts and genes (expected)
- Cells with high mt% and low total counts (dead/dying)
- Cells with very high UMIs and genes (potential doublets)
- Outlier clusters

### Plot 3: Ranked Barcode Plot
```python
# Knee plot showing cell calling
# (Usually from CellRanger or equivalent)
```

**Look for:**
- Clear inflection point
- Separation between cells and empty droplets
- Whether your filtered cells are above the knee

---

## Step 4: Apply Threshold Selection Method

Choose one method based on your data characteristics:

### Method A: Fixed Thresholds (Simple, Less Recommended)

**Use when:**
- You have strong prior knowledge of expected ranges
- You're analyzing very standard samples (e.g., PBMC 10x v3)
- You have multiple similar samples to compare

**Example thresholds:**
```python
# Conservative (strict)
min_genes = 500
max_genes = 6000
max_mt = 10

# Moderate (balanced)
min_genes = 200
max_genes = 7000
max_mt = 15

# Permissive (loose - for rare cells)
min_genes = 100
max_genes = 10000
max_mt = 20
```

**⚠️ Warning:** These are examples only. Don't copy-paste without justification!

### Method B: MAD-Based Thresholds (Recommended)

**Use when:**
- You want adaptive, data-driven thresholds
- You don't have strong priors
- You want to identify outliers objectively

**Median Absolute Deviation (MAD) Method:**

```python
import numpy as np

def is_outlier(adata, metric, n_mads=5):
    """
    Identify outliers using MAD method
    n_mads: number of MADs from median (typically 3-5)
    """
    M = adata.obs[metric]
    median = np.median(M)
    mad = np.median(np.abs(M - median))
    
    lower = median - n_mads * mad
    upper = median + n_mads * mad
    
    return (M < lower) | (M > upper)

# Apply
adata.obs['outlier_counts'] = is_outlier(adata, 'log10_total_counts', n_mads=5)
adata.obs['outlier_genes'] = is_outlier(adata, 'log10_n_genes', n_mads=5)
adata.obs['outlier_mt'] = is_outlier(adata, 'pct_counts_mt', n_mads=3)

# Filter
adata = adata[~(adata.obs.outlier_counts | 
                adata.obs.outlier_genes | 
                adata.obs.outlier_mt)]
```

**Adjust n_mads:**
- **Strict:** n_mads = 3 (removes more cells)
- **Moderate:** n_mads = 4
- **Permissive:** n_mads = 5 (keeps more cells)

### Method C: Quantile-Based (Alternative)

**Use when:**
- You want to remove a fixed percentage of extreme cells
- You have very large datasets
- MAD is too sensitive to small outlier groups

```python
# Remove bottom 1% and top 1% by total counts
lower_bound = adata.obs['total_counts'].quantile(0.01)
upper_bound = adata.obs['total_counts'].quantile(0.99)

# For mt%, typically only upper threshold
mt_upper = adata.obs['pct_counts_mt'].quantile(0.95)

adata = adata[(adata.obs['total_counts'] > lower_bound) &
              (adata.obs['total_counts'] < upper_bound) &
              (adata.obs['pct_counts_mt'] < mt_upper)]
```

---

## Step 5: Apply Tissue-Specific Guidelines

### Brain/Neuronal Tissue
```python
# Neurons have higher mt% due to high energy demands
min_genes = 200
max_mt = 20  # More permissive!
min_counts = 500
```

**Rationale:** Neurons are metabolically active; 15-20% mt% is biologically normal, not dead cells.

### PBMC / Blood
```python
# Blood cells should have low mt%
min_genes = 200
max_mt = 10  # Stricter!
min_counts = 500
```

**Rationale:** High mt% in blood usually indicates cell death during processing.

### Tumor Tissue
```python
# High heterogeneity expected
min_genes = 100  # More permissive for low-expressing cells
max_mt = 20
min_counts = 200
```

**Rationale:** Tumor cells can have diverse states; don't over-filter.

### Sorted/Enriched Populations
```python
# Expect more uniform QC distributions
# Use MAD method with lower n_mads
min_genes = 300
max_mt = 10
min_counts = 500
```

**Rationale:** Sorting enriches for specific cell types; distributions should be tighter.

### Frozen vs Fresh Tissue

**Frozen tissue typically has:**
- Higher mt% (freezing stress)
- Lower gene counts
- More ambient RNA

**Adjust thresholds accordingly:**
```python
# Frozen tissue - more permissive
max_mt = 20  # vs 10 for fresh
min_genes = 150  # vs 200 for fresh
```

---

## Step 6: Check for Over-Filtering

After applying thresholds, **validate you haven't removed biology:**

### Validation Checklist

- [ ] **Cell count reasonable?**
  - Lost < 30% of cells (unless very poor quality)
  - If > 50% filtered → reconsider thresholds
  
- [ ] **Known cell types present?**
  - Check for canonical markers of expected cell types
  - Plot marker genes on UMAP (even pre-clustering)
  
- [ ] **QC metric separation from biology?**
  - After filtering, do clusters still separate by QC metrics?
  - If yes → may need gentler filtering
  
- [ ] **Compare to published data:**
  - Similar tissue type should have similar QC distributions
  - Major discrepancies need investigation

### Red Flags (You May Have Over-Filtered)

🚩 **Lost expected rare population**
- Solution: Lower thresholds, especially for min_genes

🚩 **Removed an entire sample**
- Solution: Apply per-sample thresholds, not global

🚩 **Remaining cells have high batch/QC-driven clustering**
- Solution: Filtering didn't address underlying quality issues

🚩 **Lost > 50% of cells with moderate thresholds**
- Solution: Investigate upstream issues (poor dissociation, cell death)

---

## Step 7: Document Your Decisions

**Create a QC report documenting:**

### Required Documentation

1. **Threshold Values Used**
```
Minimum genes: 200
Maximum genes: 6000
Maximum mt%: 15
Method: MAD with n_mads=5
```

2. **Rationale**
```
Thresholds chosen based on:
- Tissue type: PBMC (expect low mt%)
- Fresh tissue (not frozen)
- MAD method to adapt to sample-specific distributions
- Visual inspection of QC plots (no bimodality observed)
```

3. **Before/After Statistics**
```
Before filtering: 5,234 cells
After filtering: 4,521 cells (13.6% removed)

Median genes before: 1,200
Median genes after: 1,450

Median mt% before: 8.5%
Median mt% after: 6.2%
```

4. **Known Marker Check**
```
Confirmed presence of:
- CD3E (T cells) ✓
- CD79A (B cells) ✓
- LYZ (Monocytes) ✓
- NKG7 (NK cells) ✓
```

5. **Figures to Include**
- Violin plots before and after filtering
- Scatter plots showing filtered cells (in different color)
- Knee plot with cell calling threshold

---

## Common Scenarios & Solutions

### Scenario 1: Bimodal n_genes Distribution

**Observation:** Two clear peaks in gene count distribution

**Possible causes:**
1. Two cell types with very different transcriptional complexity
2. Mixture of singlets and doublets
3. Empty droplets not properly filtered

**Solution:**
```python
# Don't use global min_genes threshold
# Instead, remove only extreme low-count outliers
min_genes = 100  # Very permissive

# Rely on doublet detection for high-count outliers
# Use clustering to separate cell types
```

### Scenario 2: Gradual mt% Increase, No Clear Cutoff

**Observation:** Mt% gradually increases, no clear "dead cell" population

**Possible causes:**
1. Some cell stress during dissociation
2. Tissue-specific variation
3. Biological signal (e.g., metabolic states)

**Solution:**
```python
# Use quantile approach
mt_cutoff = adata.obs['pct_counts_mt'].quantile(0.95)

# Or use MAD but make it permissive
# Pair with downstream QC (check if mt% drives clustering)
```

### Scenario 3: Sample Has Very Low Counts

**Observation:** One sample has much lower UMI counts than others

**Possible causes:**
1. Sequencing depth variation
2. Library prep failure
3. Biological difference (e.g., cell death in disease sample)

**Solution:**
```python
# Apply per-sample thresholds
for sample in adata.obs['sample'].unique():
    sample_cells = adata.obs['sample'] == sample
    
    # Calculate sample-specific MAD thresholds
    # Then apply filtering
```

**Or consider:** Exclude the entire sample if quality is too poor

### Scenario 4: High-Quality Data, Few Outliers

**Observation:** Clean, tight distributions for all QC metrics

**Solution:**
```python
# Use gentle filtering (don't over-filter good data)
# Remove only extreme outliers
# Focus on doublet detection
```

---

## QC Checklist Summary

Use this checklist for every dataset:

### Before Filtering
- [ ] Understand tissue type and biology
- [ ] Generate all QC metrics
- [ ] Create violin and scatter plots
- [ ] Check for sample-to-sample variation
- [ ] Identify distribution patterns (unimodal, bimodal, outliers)

### Choosing Thresholds
- [ ] Select appropriate method (MAD, fixed, quantile)
- [ ] Apply tissue-specific adjustments
- [ ] Consider fresh vs frozen
- [ ] Account for rare populations if expected

### After Filtering
- [ ] Check cells removed percentage (< 30% ideal)
- [ ] Validate known markers present
- [ ] Verify QC metrics don't drive clustering
- [ ] Compare to similar published datasets
- [ ] Document all decisions with rationale

### Final Validation
- [ ] Run doublet detection
- [ ] Check if further filtering needed
- [ ] Save both raw and filtered objects
- [ ] Create before/after comparison report

---

## Templates

### QC Report Template

```markdown
# QC Report: [Dataset Name]

## Sample Information
- Tissue type: ___________
- Protocol: ___________
- Number of samples: ___
- Expected cell types: ___________

## QC Metrics Summary

### Before Filtering
| Metric | Median | MAD | Min | Max |
|--------|--------|-----|-----|-----|
| n_genes | | | | |
| total_counts | | | | |
| pct_counts_mt | | | | |

Cells: _____

### Threshold Selection
- Method: [MAD / Fixed / Quantile]
- Minimum genes: ___
- Maximum genes: ___
- Maximum mt%: ___
- Rationale: ___________

### After Filtering
| Metric | Median | MAD | Min | Max |
|--------|--------|-----|-----|-----|
| n_genes | | | | |
| total_counts | | | | |
| pct_counts_mt | | | | |

Cells: _____ (___% removed)

## Validation
- [ ] Known markers detected: [list]
- [ ] No QC-driven clustering
- [ ] Compared to [reference dataset]
- [ ] All samples retained

## Figures
1. QC violin plots (before/after)
2. Scatter plots with thresholds marked
3. Per-sample QC distributions

## Conclusion
[Brief statement on data quality and suitability for downstream analysis]
```

---

## Key Takeaways

1. **No universal thresholds** - Always consider biological context
2. **Visualize first** - Plots reveal distribution patterns
3. **MAD method preferred** - Adaptive and data-driven
4. **Tissue-specific rules** - Brain ≠ Blood ≠ Tumor
5. **Validate don't remove biology** - Check for known markers
6. **Document everything** - Reproducibility requires justification
7. **Iterate if needed** - QC is not one-and-done
8. **Save raw data** - Never overwrite; filter creates new object

**Remember:** The goal is to enable good downstream analysis, not to achieve perfect-looking QC plots at the expense of biological signal.

---

**Related Resources:**
- `../course_single_cell_processing.md` - Full course context
- `../labs/lab08_qc_metrics.ipynb` - Hands-on QC practice
- `../labs/lab09_filtering_doublets.ipynb` - Filtering implementation
- Recommended reading: Luecken & Theis (2019) "Current best practices in single-cell RNA-seq analysis"
