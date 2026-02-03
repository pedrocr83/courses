# Module 11: Trajectory Validation Strategies

**How to validate and gain confidence in inferred trajectories**

---

## Overview

Trajectory inference is **computational prediction**, not ground truth. This module teaches you systematic validation strategies to:
1. Assess reliability of inferred trajectories
2. Identify artifacts vs real biological ordering
3. Build evidence for (or against) trajectory claims
4. Appropriately report trajectory results with caveats

**Duration:** 2 hours  
**Prerequisites:** Modules 2-6 (trajectory inference methods)

---

## Learning Objectives

By the end of this module, you will be able to:
- [ ] Design multi-pronged validation strategies for trajectories
- [ ] Use marker genes to validate pseudotime ordering
- [ ] Assess sensitivity of trajectories to method and parameter choices
- [ ] Identify common trajectory inference artifacts
- [ ] Report trajectory results with appropriate confidence levels

---

## Why Validation Matters

### The Problem
Trajectory inference tools **will always** produce a trajectory, even when:
- Cells are actually discrete types (no true continuity)
- Batch effects drive the ordering (not biology)
- Clusters are arbitrarily connected
- Multiple valid orderings exist

### The Solution
**Multi-level validation** using orthogonal evidence:
1. **Internal consistency**: Method agreement, parameter robustness
2. **Marker genes**: Known temporal genes change appropriately
3. **External data**: Comparison with time-course, lineage tracing
4. **Perturbation experiments**: Predicted drivers validated functionally

---

## Validation Strategy Framework

### Level 1: Internal Validation (Required)

**Can be done with your current data alone**

#### 1.1 Marker Gene Validation
Check if known early/late markers change in expected direction

#### 1.2 Method Agreement
Compare results across 2-3 trajectory inference methods

#### 1.3 Parameter Sensitivity
Test if conclusions hold across reasonable parameter ranges

#### 1.4 Biological Plausibility
Assess if trajectory makes biological sense

---

### Level 2: Cross-Dataset Validation (Recommended)

**Requires related datasets**

#### 2.1 Published Trajectory Comparison
Compare to published trajectories in same biological system

#### 2.2 Time-Course Integration
Align pseudotime with real time-course data (if available)

#### 2.3 Cross-Species Conservation
Check if trajectory structure is conserved across species

---

### Level 3: Experimental Validation (Gold Standard)

**Requires follow-up experiments**

#### 3.1 Lineage Tracing
Direct observation of cell fate decisions

#### 3.2 Perturbation Validation
Knockdown/overexpression of predicted drivers

#### 3.3 Temporal Sampling
Time-series experiment matching pseudotime stages

#### 3.4 Protein-Level Validation
Flow cytometry, immunofluorescence of key markers

---

## Detailed Validation Methods

### Method 1: Marker Gene Validation

**Goal:** Verify known temporal genes change in expected direction

**Implementation:**

```python
import scanpy as sc
import numpy as np

# 1. Define known early/late markers (from literature)
early_markers = ['NANOG', 'POU5F1', 'SOX2']  # Example: stem cell
late_markers = ['EOMES', 'T', 'MESP1']       # Example: mesoderm

# 2. Calculate correlation of markers with pseudotime
from scipy.stats import spearmanr

results = []
for gene in early_markers + late_markers:
    if gene in adata.var_names:
        gene_expr = adata[:, gene].X.toarray().flatten()
        corr, pval = spearmanr(gene_expr, adata.obs['dpt_pseudotime'])
        results.append({
            'gene': gene,
            'type': 'early' if gene in early_markers else 'late',
            'correlation': corr,
            'p_value': pval
        })

# 3. Expected pattern:
# Early markers: negative correlation (high at start, low at end)
# Late markers: positive correlation (low at start, high at end)

import pandas as pd
df = pd.DataFrame(results)
print(df)

# 4. Visualize
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

for markers, ax, title in zip([early_markers, late_markers], axes, ['Early Markers', 'Late Markers']):
    for gene in markers:
        if gene in adata.var_names:
            gene_expr = adata[:, gene].X.toarray().flatten()
            ax.scatter(adata.obs['dpt_pseudotime'], gene_expr, alpha=0.3, s=5, label=gene)
    ax.set_xlabel('Pseudotime')
    ax.set_ylabel('Expression')
    ax.set_title(title)
    ax.legend()
plt.tight_layout()
plt.show()
```

**Interpretation:**
- ✅ **PASS:** Early markers decrease, late markers increase along pseudotime
- ⚠️ **CAUTION:** Some markers behave as expected, others don't
- ❌ **FAIL:** No clear pattern or opposite of expected

---

### Method 2: Root Cell Validation

**Goal:** Verify root cell selection is biologically appropriate

**Implementation:**

```python
# 1. Identify root cell(s)
root_cell = adata.obs['dpt_pseudotime'].idxmin()

# 2. Check marker expression in root
root_markers = adata[root_cell, early_markers].X.toarray()
print(f"Root cell expression of early markers:\n{root_markers}")

# 3. Test multiple plausible roots
import numpy as np

# Find cells in earliest cluster
early_cluster = ...  # Identify based on markers
early_cells = adata[adata.obs['leiden'] == early_cluster].obs_names

# Test each as potential root
root_tests = []
for test_root in np.random.choice(early_cells, size=min(10, len(early_cells)), replace=False):
    # Re-run pseudotime with new root
    adata_test = adata.copy()
    adata_test.uns['iroot'] = np.where(adata_test.obs_names == test_root)[0][0]
    sc.tl.dpt(adata_test)
    
    # Calculate marker correlation
    early_corr = [spearmanr(adata_test[:, g].X.toarray().flatten(), 
                           adata_test.obs['dpt_pseudotime'])[0] 
                  for g in early_markers if g in adata.var_names]
    
    root_tests.append({
        'root': test_root,
        'mean_early_corr': np.mean(early_corr)
    })

# 4. Check sensitivity to root choice
import pandas as pd
df = pd.DataFrame(root_tests)
print(f"Root sensitivity:\n{df.describe()}")

# Trajectory is robust if correlations are similar across roots
```

**Interpretation:**
- ✅ **ROBUST:** Results consistent across plausible root choices
- ⚠️ **SENSITIVE:** Results change moderately with root
- ❌ **UNSTABLE:** Results highly dependent on root (trajectory unreliable)

---

### Method 3: Method Agreement

**Goal:** Check if multiple methods agree on trajectory structure

**Implementation:**

```python
# Run 3 different methods on same data

# Method 1: Diffusion pseudotime
sc.pp.neighbors(adata, n_pcs=30)
sc.tl.diffmap(adata)
sc.tl.dpt(adata)
adata.obs['dpt_pt'] = adata.obs['dpt_pseudotime']

# Method 2: PAGA + diffusion
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata)
# (diffusion pseudotime already calculated)

# Method 3: Palantir or Monocle (if available)
# import palantir
# pr_res = palantir.core.run_palantir(adata, ...)
# adata.obs['palantir_pt'] = pr_res.pseudotime

# Compare pseudotime orderings
from scipy.stats import spearmanr

corr_dpt_paga = spearmanr(adata.obs['dpt_pt'], adata.obs['dpt_pseudotime'])[0]
print(f"DPT vs PAGA pseudotime correlation: {corr_dpt_paga:.3f}")

# Visualize agreement
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(adata.obs['dpt_pt'], adata.obs['palantir_pt'], 
           alpha=0.3, s=5)
ax.set_xlabel('DPT Pseudotime')
ax.set_ylabel('Palantir Pseudotime')
ax.set_title(f"Method Agreement (r={corr_dpt_paga:.3f})")
plt.show()

# Check if cell ordering is similar
from scipy.stats import kendalltau
tau, pval = kendalltau(adata.obs['dpt_pt'], adata.obs['palantir_pt'])
print(f"Kendall's tau: {tau:.3f} (p={pval:.2e})")
```

**Interpretation:**
- ✅ **STRONG AGREEMENT:** r > 0.8, tau > 0.7
- ⚠️ **MODERATE AGREEMENT:** r = 0.5-0.8, tau = 0.4-0.7
- ❌ **POOR AGREEMENT:** r < 0.5 (methods disagree fundamentally)

---

### Method 4: Parameter Sensitivity Analysis

**Goal:** Ensure conclusions aren't driven by arbitrary parameter choices

**Implementation:**

```python
# Test multiple parameter combinations
import pandas as pd
from itertools import product

# Parameters to test
n_neighbors_values = [10, 20, 30, 50]
n_pcs_values = [20, 30, 40]

results = []
for n_neighbors, n_pcs in product(n_neighbors_values, n_pcs_values):
    # Re-run analysis
    adata_test = adata.copy()
    sc.pp.neighbors(adata_test, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.diffmap(adata_test)
    sc.tl.dpt(adata_test)
    
    # Calculate marker correlations
    early_corr = np.mean([spearmanr(adata_test[:, g].X.toarray().flatten(), 
                                    adata_test.obs['dpt_pseudotime'])[0]
                         for g in early_markers if g in adata.var_names])
    
    late_corr = np.mean([spearmanr(adata_test[:, g].X.toarray().flatten(), 
                                   adata_test.obs['dpt_pseudotime'])[0]
                        for g in late_markers if g in adata.var_names])
    
    results.append({
        'n_neighbors': n_neighbors,
        'n_pcs': n_pcs,
        'early_marker_corr': early_corr,
        'late_marker_corr': late_corr,
        'consistent': (early_corr < -0.3 and late_corr > 0.3)  # Expected pattern
    })

# Summary
df = pd.DataFrame(results)
print(df)
print(f"\nConsistent across {df['consistent'].sum()} / {len(df)} parameter combinations")

# Visualize
import seaborn as sns
pivot = df.pivot(index='n_neighbors', columns='n_pcs', values='late_marker_corr')
sns.heatmap(pivot, annot=True, cmap='RdBu_r', center=0, vmin=-1, vmax=1)
plt.title('Late Marker Correlation Across Parameters')
plt.show()
```

**Interpretation:**
- ✅ **ROBUST:** Consistent pattern across >80% of parameter combinations
- ⚠️ **MODERATELY ROBUST:** Consistent in 50-80% of combinations
- ❌ **FRAGILE:** Pattern changes frequently with parameters

---

### Method 5: Batch Effect Check

**Goal:** Ensure trajectory isn't driven by batch/technical effects

**Implementation:**

```python
# Check if pseudotime correlates with batch
if 'batch' in adata.obs.columns:
    # Visualize
    sc.pl.umap(adata, color=['dpt_pseudotime', 'batch'], ncols=2)
    
    # Statistical test
    from scipy.stats import kruskal
    
    batches = adata.obs['batch'].unique()
    pt_by_batch = [adata[adata.obs['batch'] == b].obs['dpt_pseudotime'].values 
                   for b in batches]
    
    stat, pval = kruskal(*pt_by_batch)
    print(f"Kruskal-Wallis test for pseudotime vs batch: p={pval:.2e}")
    
    if pval < 0.05:
        print("⚠️ WARNING: Pseudotime significantly differs by batch!")
        print("Trajectory may be driven by batch effects, not biology.")
    
    # Check if early/late stages are batch-specific
    early = adata[adata.obs['dpt_pseudotime'] < np.percentile(adata.obs['dpt_pseudotime'], 25)]
    late = adata[adata.obs['dpt_pseudotime'] > np.percentile(adata.obs['dpt_pseudotime'], 75)]
    
    print("\nBatch composition:")
    print("Early stage:", early.obs['batch'].value_counts(normalize=True))
    print("Late stage:", late.obs['batch'].value_counts(normalize=True))
```

**Interpretation:**
- ✅ **CLEAN:** No batch-pseudotime association (p > 0.05)
- ⚠️ **CONFOUNDED:** Moderate batch effect (re-integrate before trajectory)
- ❌ **BATCH-DRIVEN:** Strong batch effect (trajectory is artifact)

---

### Method 6: Biological Plausibility Assessment

**Checklist for assessing biological sensibility:**

- [ ] **Trajectory makes biological sense**
  - Is the proposed process known to be continuous?
  - Are the cell types along the trajectory related?
  - Does the directionality match known biology?

- [ ] **Cell types are appropriately ordered**
  - Progenitors before differentiated cells?
  - Intermediate states exist and make sense?
  - Terminal states are biologically terminal?

- [ ] **Gene expression patterns are coherent**
  - Genes with known temporal roles behave correctly?
  - Gene modules change gradually (not abruptly)?
  - Pathway activation follows expected sequence?

- [ ] **Trajectory structure is supported**
  - Branching points have biological basis?
  - Number of branches matches known fates?
  - Rare intermediates detectable (if expected)?

**Red Flags (Biologically Implausible):**
- 🚩 Mature cells "before" progenitors
- 🚩 Unrelated cell types on same trajectory
- 🚩 Arbitrary cluster ordering
- 🚩 Trajectory contradicts known lineage relationships

---

## Validation Report Template

### Trajectory Validation Checklist

**Dataset:** ______________________  
**Biological Process:** ______________________  
**Methods Used:** ______________________

---

#### 1. Marker Gene Validation

| Marker | Type | Expected Direction | Observed Correlation | PASS/FAIL |
|--------|------|-------------------|----------------------|-----------|
| Gene1  | Early| Negative          | -0.72               | ✅ PASS    |
| Gene2  | Late | Positive          | +0.81               | ✅ PASS    |
| ...    | ...  | ...               | ...                 | ...       |

**Summary:** ___/___markers behave as expected

---

#### 2. Root Cell Validation

**Root cluster:** ___________  
**Root marker expression:** [High/Low] for expected markers  
**Root sensitivity:** [Robust/Sensitive/Unstable]  
**Conclusion:** ✅ / ⚠️ / ❌

---

#### 3. Method Agreement

| Method Pair | Spearman r | Kendall τ | Agreement |
|-------------|------------|-----------|-----------|
| DPT vs PAGA | 0.85       | 0.73      | ✅ Strong  |
| DPT vs Monocle| 0.72     | 0.61      | ⚠️ Moderate|
| ...         | ...        | ...       | ...       |

**Conclusion:** [Strong/Moderate/Poor] agreement across methods

---

#### 4. Parameter Sensitivity

**Parameters tested:** _____________  
**Consistent results:** ___% of combinations  
**Conclusion:** [Robust/Moderate/Fragile]

---

#### 5. Batch Effect Check

**Batch-pseudotime association:** p = _____  
**Batch composition balance:** [Balanced/Imbalanced]  
**Conclusion:** [Clean/Confounded/Batch-driven]

---

#### 6. Biological Plausibility

- [ ] Process known to be continuous: YES / NO
- [ ] Cell type ordering makes sense: YES / NO
- [ ] Marker patterns coherent: YES / NO
- [ ] Trajectory structure supported: YES / NO

**Red flags identified:** ________________

---

### Overall Confidence Assessment

**Confidence Level:** [HIGH/MEDIUM/LOW]

**Rationale:** 
_______________________________________

**Caveats:**
_______________________________________

**Recommended Follow-Up Validation:**
1. _______________________________________
2. _______________________________________
3. _______________________________________

---

## Common Validation Failures & Solutions

### Failure 1: Markers Don't Follow Expected Pattern

**Symptoms:**
- Early markers don't decrease along pseudotime
- Late markers don't increase
- Mixed or no clear pattern

**Possible Causes:**
1. Wrong root selection
2. Trajectory doesn't actually exist (discrete types)
3. Wrong markers (check literature)
4. Batch effects driving ordering

**Solutions:**
```python
# Test different roots
# Check for discrete types (try clustering instead)
# Re-check marker literature
# Correct batch effects before trajectory inference
```

---

### Failure 2: Methods Strongly Disagree

**Symptoms:**
- Low correlation between pseudotime from different methods
- Contradictory ordering
- Different trajectory structures

**Possible Causes:**
1. Weak trajectory signal (multiple valid orderings)
2. Parameter-dependent results
3. Method-specific biases

**Solutions:**
```python
# Focus on robust features (present in all methods)
# Report consensus only
# Consider that trajectory may not be well-defined
# Lower confidence in conclusions
```

---

### Failure 3: Batch-Driven Trajectory

**Symptoms:**
- Pseudotime correlates with batch
- Early/late stages enriched for specific batches
- Batch separation on UMAP follows pseudotime gradient

**Solutions:**
```python
# Integrate batches BEFORE trajectory inference
import scanpy as sc
import harmonypy as hm

# Batch correction
hm_res = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
adata.obsm['X_pca_harmony'] = hm_res.Z_corr.T

# Re-build neighbors on corrected PCA
sc.pp.neighbors(adata, use_rep='X_pca_harmony')

# Re-run trajectory
sc.tl.diffmap(adata)
sc.tl.dpt(adata)

# Re-validate
```

---

## External Validation Strategies

### Strategy 1: Time-Course Integration

**If you have time-series data:**

```python
# Align pseudotime with real time
import numpy as np
from scipy.stats import spearmanr

# Add real time to metadata
adata.obs['real_time'] = ...  # In hours/days

# Calculate correlation
corr, pval = spearmanr(adata.obs['dpt_pseudotime'], adata.obs['real_time'])
print(f"Pseudotime vs Real Time: r={corr:.3f}, p={pval:.2e}")

# Expected: Positive correlation
# Strong correlation (r>0.7) = good validation
```

---

### Strategy 2: Published Trajectory Comparison

**Compare to literature:**

```python
# Load published marker list
published_early = ['GENE1', 'GENE2', ...]
published_late = ['GENE3', 'GENE4', ...]

# Check agreement
early_corr = [spearmanr(adata[:, g].X.toarray().flatten(), 
                       adata.obs['dpt_pseudotime'])[0]
              for g in published_early if g in adata.var_names]

print(f"Agreement with published early markers: {np.mean(early_corr):.3f}")
# Expected: Negative (early markers decrease)

# Can also visually compare trajectory structures
```

---

### Strategy 3: Experimental Validation Design

**Propose follow-up experiments:**

1. **Lineage Tracing**
   - Use genetic barcoding or viral labeling
   - Track cell fates over time
   - Confirm predicted lineage relationships

2. **Perturbation Screen**
   - Identify genes changing along pseudotime
   - Knockdown/overexpress candidate drivers
   - Observe impact on differentiation

3. **Temporal Sampling**
   - Design experiment matching pseudotime stages
   - Collect samples at predicted early/intermediate/late points
   - Validate marker expression with protein-level assays

4. **Single-Cell Resolution Validation**
   - FISH or immunofluorescence of key markers
   - Confirm spatial organization matches trajectory
   - Validate rare intermediate populations

---

## Reporting Trajectories in Publications

### What to Include

1. **Methods Section**
   - Trajectory inference method(s) used
   - Parameter values
   - Root selection criteria
   - Software versions

2. **Results Section**
   - Trajectory structure (with figure)
   - Marker validation results
   - Method agreement (if applicable)
   - Confidence level statement

3. **Supplementary**
   - Full validation checklist
   - Parameter sensitivity analysis
   - Alternative methods comparison
   - Marker gene trends along pseudotime

---

### Example Statements

**HIGH CONFIDENCE:**
> "We inferred developmental trajectory using diffusion pseudotime (DPT). The trajectory was validated by: (1) expected expression patterns of known early (SOX2, NANOG) and late (EOMES, T) markers (Supplementary Fig. X), (2) agreement with an independent method (Monocle 3, r=0.85), and (3) consistency across parameter ranges (n_neighbors=10-50). Results were robust to root cell selection (Supplementary Fig. Y)."

**MEDIUM CONFIDENCE:**
> "We inferred a potential differentiation trajectory using DPT. While several known late markers (e.g., EOMES) increased along pseudotime, early marker patterns were less consistent. The trajectory structure showed moderate agreement with PAGA (r=0.68) but was sensitive to root selection. These results suggest a general developmental progression but require experimental validation."

**LOW CONFIDENCE:**
> "Computational trajectory inference suggested a possible ordering of cell states. However, validation markers showed mixed patterns, and results varied substantially between methods. The inferred trajectory should be interpreted cautiously and requires orthogonal validation before biological interpretation."

---

## Summary: Validation Best Practices

### DO:
✅ Validate with multiple orthogonal approaches  
✅ Test robustness to methods and parameters  
✅ Use known markers extensively  
✅ Check for batch confounding  
✅ Report confidence levels honestly  
✅ Propose experimental follow-up  
✅ Show validation results in supplementary  

### DON'T:
❌ Trust a single method without validation  
❌ Ignore method disagreements  
❌ Claim high confidence without evidence  
❌ Skip marker validation  
❌ Assume trajectory exists because algorithm produced one  
❌ Report trajectories without caveats  
❌ Ignore batch effects  

---

## Practice Exercise

**Use your own trajectory analysis and complete the validation checklist:**

1. [ ] Identify 5+ known temporal markers from literature
2. [ ] Calculate correlations with pseudotime
3. [ ] Run at least 2 trajectory methods
4. [ ] Compare results (calculate correlation)
5. [ ] Test 3+ parameter combinations
6. [ ] Check for batch effects
7. [ ] Complete validation report template
8. [ ] Assign confidence level (HIGH/MEDIUM/LOW)
9. [ ] Propose 2-3 experimental validations

---

## Resources

### Marker Databases
- **CellMarker:** Developmental stage markers
- **Enrichr:** Pathway databases for stage-specific genes
- **PanglaoDB:** Cell type temporal markers

### Validation Papers
- Saelens et al. (2019) "A comparison of single-cell trajectory inference methods"
- Tritschler et al. (2019) "Concepts and limitations for learning developmental trajectories"
- Street et al. (2018) "Slingshot: cell lineage and pseudotime inference"

### Tools
- **scvelo:** RNA velocity validation
- **CellRank:** Fate probability validation
- **tradeSeq:** Statistical testing along trajectories

---

**Next Steps:**
- Complete validation checklist for your data
- Review validation section in `../assignments/assignment1_pseudotime.md`
- See examples in `../labs/lab03_diffusion_pt.ipynb`

**Remember:** A trajectory without validation is just a pretty picture. Build the evidence!
