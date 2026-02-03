# Integration Method Decision Guide

**Interactive decision tree and comprehensive guide for choosing the right batch correction method**

---

## Quick Decision Tree

Answer these questions to get a recommended method:

### Question 1: Dataset Size
**How many cells total across all batches?**

- **A) < 10,000 cells** → Go to Q2a
- **B) 10,000 - 100,000 cells** → Go to Q2b
- **C) > 100,000 cells** → Go to Q2c

---

### Q2a: Small Dataset (< 10K cells)

**Do you have GPU access?**
- **Yes** → **Recommended: scVI** (excellent performance, quick on small data)
- **No** → Go to Q3a

#### Q3a: Number of batches?
- **2-3 batches** → **Recommended: Seurat CCA** (reliable for few batches)
- **4+ batches** → **Recommended: Harmony** (handles many batches well)

---

### Q2b: Medium Dataset (10-100K cells)

**Primary consideration?**
- **A) Speed** → **Recommended: Harmony** (fastest for this size)
- **B) Accuracy** → Go to Q3b
- **C) Interpretability** → **Recommended: Seurat CCA** (understand anchors)

#### Q3b: Do you have GPU?
- **Yes** → **Recommended: scVI** (excellent accuracy)
- **No** → **Recommended: Harmony or Seurat RPCA**

---

### Q2c: Large Dataset (> 100K cells)

**Do you have GPU with ≥16GB VRAM?**
- **Yes** → **Recommended: scVI** (scales well with GPU)
- **No** → **Recommended: Harmony** (CPU-friendly, fast)

**Alternative:** Consider downsampling for exploration, then scale up with Harmony

---

## Detailed Method Comparison

### Overview Table

| Method | Speed | Accuracy | Scalability | GPU Needed | Interpretability | Best For |
|--------|-------|----------|-------------|------------|------------------|----------|
| **Harmony** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | No | ⭐⭐⭐ | Large datasets, many batches |
| **Seurat CCA** | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | No | ⭐⭐⭐⭐⭐ | Few batches, moderate size |
| **Seurat RPCA** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | No | ⭐⭐⭐⭐ | Query-to-reference mapping |
| **MNN (fastMNN)** | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | No | ⭐⭐⭐ | Limited shared cell types |
| **scVI** | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Yes* | ⭐⭐ | Complex batch structure, multi-condition |
| **scANVI** | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Yes* | ⭐⭐ | When you have some labels |

*Can run on CPU but very slow

---

## Method Profiles

### Harmony

**When to use:**
- ✅ Large datasets (100K+ cells)
- ✅ Many batches (5-20+)
- ✅ Need fast results
- ✅ Limited computational resources
- ✅ Multi-condition experiments

**When NOT to use:**
- ❌ Only 2 batches (overkill, use simpler methods)
- ❌ Need reference-query mapping (use RPCA instead)
- ❌ Batch perfectly confounded with biology (nothing will save you)

**Key parameters to tune:**
```python
# Python (harmonypy)
harmony_out = harmony.run_harmony(
    data,
    meta_data,
    vars_use=['batch'],
    theta=[2],  # Higher = stronger correction (default 2)
    lambda_val=[1],  # Ridge regression parameter
    nclust=50,  # Number of clusters (default)
    max_iter=10  # More for stubborn batches
)
```

**Pros:**
- Very fast
- Scales to millions of cells
- Handles many batches effortlessly
- Preserves cell type separation well
- Works in PCA space (interpretable)

**Cons:**
- Can over-correct if theta too high
- Less control than anchor-based methods
- Black box optimization

**Typical runtime:** 1-5 minutes for 50K cells on laptop

---

### Seurat CCA (Canonical Correlation Analysis)

**When to use:**
- ✅ 2-5 batches
- ✅ 10K-100K cells
- ✅ Want to understand integration (anchors are interpretable)
- ✅ Batches have substantial differences
- ✅ Need to preserve subtle cell type differences

**When NOT to use:**
- ❌ >10 batches (slow)
- ❌ Very large datasets (>200K)
- ❌ Batches are already quite similar (use faster method)

**Key parameters:**
```r
# R (Seurat v5)
integrated <- IntegrateLayers(
  object = seurat_obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  dims = 1:30,  # Number of PCs to use
  k.anchor = 5,  # Anchors per cell (higher = more correction)
  k.filter = 200,  # Min anchors per cell
  k.weight = 100  # Weight neighbors
)
```

**Pros:**
- Well-established and trusted
- Interpretable (can inspect anchors)
- Good preservation of cell type separation
- Flexible with multi-condition designs

**Cons:**
- Slow for large datasets
- Memory intensive
- Requires careful parameter tuning
- Can over-correct if not careful

**Typical runtime:** 10-30 minutes for 50K cells

---

### Seurat RPCA (Reciprocal PCA)

**When to use:**
- ✅ Reference-query mapping scenarios
- ✅ One batch is trusted reference, others are queries
- ✅ Moderate to large datasets
- ✅ Faster than CCA but similar logic

**When NOT to use:**
- ❌ All batches are equally important (use CCA or Harmony)
- ❌ No clear reference dataset

**Key parameters:**
```r
# R (Seurat v5)
integrated <- IntegrateLayers(
  object = seurat_obj,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  k.anchor = 5,
  reference = 1  # Specify reference dataset index
)
```

**Pros:**
- Faster than CCA
- Good for asymmetric integration
- Preserves reference structure
- Scales better than CCA

**Cons:**
- Requires choosing a reference
- Less flexible than CCA
- Query batches conform to reference (may lose signal)

**Typical runtime:** 5-20 minutes for 50K cells

---

### MNN / fastMNN (Mutual Nearest Neighbors)

**When to use:**
- ✅ Batches have limited shared cell types
- ✅ Want local correction (preserve global structure)
- ✅ Complex batch effects
- ✅ Moderate datasets (10K-50K)

**When NOT to use:**
- ❌ No shared cell types between batches
- ❌ Very large datasets (slow)
- ❌ Need speed

**Key parameters:**
```r
# R (batchelor package)
library(batchelor)
corrected <- fastMNN(
  sce_list,
  k = 20,  # Number of nearest neighbors
  d = 50,  # Number of PCs
  cos.norm = TRUE,  # Cosine normalization
  subset.row = hvgs  # Use HVGs only
)
```

**Pros:**
- Preserves biological heterogeneity well
- Good for complex batch effects
- Published benchmarks show good performance
- Local correction (less global over-correction)

**Cons:**
- Slow for large data
- Requires shared populations (fails if none)
- Can create artificial trajectories between batches
- Memory intensive

**Typical runtime:** 10-40 minutes for 50K cells

---

### scVI (Single-Cell Variational Inference)

**When to use:**
- ✅ Complex batch structure (batch + donor + technology)
- ✅ Have GPU (CUDA-enabled)
- ✅ Want state-of-the-art accuracy
- ✅ Multi-condition experiments
- ✅ Willing to invest time in setup

**When NOT to use:**
- ❌ No GPU (painfully slow on CPU)
- ❌ Need interpretability (black box)
- ❌ Very simple batch structure (overkill)
- ❌ Tight deadline (training takes time)

**Key parameters:**
```python
# Python (scvi-tools)
import scvi

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    continuous_covariate_keys=["pct_counts_mt"],  # Optional
)

model = scvi.model.SCVI(
    adata,
    n_layers=2,  # Neural network depth
    n_latent=30,  # Latent space dimensions
    gene_likelihood="nb"  # Negative binomial
)

model.train(
    max_epochs=400,  # Training iterations
    early_stopping=True,
    early_stopping_patience=20
)

adata.obsm["X_scvi"] = model.get_latent_representation()
```

**Pros:**
- Best accuracy in benchmarks
- Handles complex batch structures
- Can integrate counts directly (no normalization needed)
- Principled statistical model
- Supports covariates (continuous batch effects)

**Cons:**
- Requires GPU for reasonable speed
- Black box (hard to interpret)
- Longer runtime (training + hyperparameter tuning)
- Steeper learning curve
- Results can vary between runs

**Typical runtime:** 20-60 minutes for 50K cells (GPU), hours on CPU

---

### scANVI (Semi-supervised scVI)

**When to use:**
- ✅ Have cell type labels for some cells
- ✅ Want to transfer labels to unlabeled cells
- ✅ Complex batch + annotation task
- ✅ Have GPU

**When NOT to use:**
- ❌ No labels available (use scVI instead)
- ❌ Labels are unreliable
- ❌ No GPU

**Key parameters:**
```python
# Python (scvi-tools)
import scvi

# After scVI setup...
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata,
    labels_key="cell_type",
    unlabeled_category="Unknown"
)

scanvi_model.train(max_epochs=200)
```

**Pros:**
- All benefits of scVI
- Leverages known labels
- Can predict labels for new cells
- Guides integration with biological knowledge

**Cons:**
- All limitations of scVI
- Requires good quality labels for training
- Even more complex than scVI

---

## Special Scenarios

### Scenario 1: Batch Confounded with Condition

**Problem:** Batch perfectly aligns with biological condition (e.g., all control samples in batch1, all treated in batch2)

**Solutions:**
1. **DON'T integrate!** - Can't distinguish technical from biological
2. **Accept limitations** - Report results with strong caveats
3. **Redesign experiment** - If possible, re-sequence with balanced design
4. **Analyze separately** - Compare per-batch, look for consistent patterns

**What NOT to do:**
❌ Integrate anyway and claim biological differences
❌ Correct for batch and assume remainder is biological
❌ Ignore the confounding

---

### Scenario 2: Query-to-Reference Mapping

**Problem:** Have a trusted reference atlas, want to map new query data

**Recommended:** **Seurat RPCA** or **scANVI**

```r
# Seurat RPCA with reference
integrated <- IntegrateLayers(
  object = combined_obj,
  method = RPCAIntegration,
  reference = 1,  # Reference dataset index
  dims = 1:30
)
```

**Why:** Preserves reference structure, efficiently maps queries

---

### Scenario 3: Cross-Technology Integration

**Problem:** Integrating 10x with Smart-seq2, or human with mouse

**Recommended:** **scVI** or **Harmony**

**Special considerations:**
- Cross-technology: Different gene coverage → use shared genes only
- Cross-species: Need homology mapping first
- May require stronger correction (increase theta/epochs)

```python
# scVI with strong batch correction
model = scvi.model.SCVI(adata, n_layers=3)  # Deeper network
model.train(max_epochs=600)  # More training
```

---

### Scenario 4: Very Simple Batch Effect

**Problem:** 2 batches, very similar, minor technical variation

**Recommended:** **Linear regression correction** or **ComBat**

**Why:** Don't over-complicate. Simple methods often sufficient.

```python
# Scanpy - regress out batch
sc.pp.regress_out(adata, ['batch'])
```

---

### Scenario 5: Batch + Donor Structure

**Problem:** Multiple batches AND multiple donors (nested structure)

**Recommended:** **scVI** (can handle multiple covariates) or **Harmony** with both covariates

```python
# Harmony with multiple variables
harmony_out = harmony.run_harmony(
    data,
    meta_data,
    vars_use=['batch', 'donor'],  # Both covariates
    theta=[2, 2]  # Separate theta for each
)

# scVI
scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="batch",
    categorical_covariate_keys=["donor"]  # Additional covariate
)
```

---

## Evaluation Checklist

After integration, evaluate using these criteria:

### Quantitative Metrics

#### Batch Mixing (Higher = Better)
- [ ] **iLISI (Inverse Simpson Index)** - Should approach # of batches
- [ ] **kBET (k-nearest Batch Effect Test)** - Acceptance rate > 0.7
- [ ] **Graph connectivity** - Batches well-mixed in nearest neighbor graph

#### Biological Conservation (Higher = Better)
- [ ] **cLISI (Cell-type LISI)** - Should stay low (cells of same type together)
- [ ] **ARI (Adjusted Rand Index)** - Cluster agreement before/after
- [ ] **NMI (Normalized Mutual Information)** - Information preservation
- [ ] **ASW (Silhouette Width)** - Cell type separation

### Visual Inspection

- [ ] **UMAP by batch:** Batches should mix within cell types
- [ ] **UMAP by cell type:** Cell types should stay separated
- [ ] **Marker genes:** Canonical markers still identify expected types
- [ ] **Cluster composition:** No batch-dominated clusters

### Biological Validation

- [ ] **Known markers preserved:** Expected cell types detectable
- [ ] **No artificial mixing:** Cell types that shouldn't mix, don't
- [ ] **Biological variation retained:** Condition differences preserved (if applicable)

---

## Parameter Tuning Guide

### If Integration Is Too Weak (Under-Correction)
**Symptoms:** Batches still clearly separated, same-cell-types don't mix

**Solutions:**
- **Harmony:** Increase `theta` (try 3-5)
- **Seurat:** Increase `k.anchor`, decrease `k.filter`
- **scVI:** Increase `n_layers` or `max_epochs`
- **General:** Use more PCs, include more HVGs

### If Integration Is Too Strong (Over-Correction)
**Symptoms:** Cell types collapse, markers lose specificity, biology lost

**Solutions:**
- **Harmony:** Decrease `theta` (try 1-1.5)
- **Seurat:** Decrease `k.anchor`, increase `k.filter`
- **scVI:** Decrease `n_latent`, stop training earlier
- **General:** Use fewer PCs, integrate within cell types first

---

## Step-by-Step Workflow

### 1. Pre-Integration
```python
# Ensure all batches processed identically
# Same normalization, same HVG selection
# Same PCA (if applicable)
```

### 2. Integration
```python
# Choose method based on decision tree
# Start with default parameters
# Save the integration object
```

### 3. Evaluation
```python
# Compute metrics (LISI, ARI, etc.)
# Generate UMAPs (by batch, by cell type)
# Check marker genes
```

### 4. Iteration (if needed)
```python
# If over/under-corrected, adjust parameters
# Compare multiple methods
# Document what you tried
```

### 5. Downstream Analysis
```python
# Use integrated embedding for:
# - Clustering
# - UMAP visualization
# - Differential expression (use original counts!)
```

---

## When NOT to Integrate

Sometimes integration is NOT the answer:

### Don't Integrate If:
1. **Batch confounded with biology** (can't disentangle)
2. **Truly different biology** (e.g., healthy vs disease, should be different!)
3. **No shared cell types** (nothing to integrate)
4. **Very high quality data** with minimal batch effects
5. **Exploratory phase** (analyze separately first to understand differences)

### Alternatives to Integration:
- **Analyze separately** and compare results
- **Meta-analysis** across datasets
- **Reference mapping** instead of full integration
- **Batch as biological question** (test for batch effects explicitly)

---

## Quick Reference Card

```
┌─────────────────────────────────────────────────────────┐
│           INTEGRATION METHOD QUICK SELECTOR             │
├─────────────────────────────────────────────────────────┤
│ NEED SPEED + MANY BATCHES      → Harmony               │
│ NEED ACCURACY + HAVE GPU        → scVI                  │
│ MODERATE SIZE + INTERPRETABLE   → Seurat CCA           │
│ QUERY-TO-REFERENCE              → Seurat RPCA          │
│ LIMITED SHARED CELL TYPES       → fastMNN              │
│ HAVE SOME LABELS                → scANVI               │
│ SIMPLE BATCH EFFECT             → Regression/ComBat    │
└─────────────────────────────────────────────────────────┘
```

---

## Recommended Reading

- Luecken et al. (2022) "Benchmarking atlas-level data integration" *Nat Methods*
- Tran et al. (2020) "A benchmark of batch-effect correction methods" *Genome Biol*
- Korsunsky et al. (2019) "Fast, sensitive and accurate integration" *Nat Methods* (Harmony)
- Lopez et al. (2018) "Deep generative modeling for single-cell transcriptomics" *Nat Methods* (scVI)

---

**Remember:** 
- No method is perfect for all datasets
- **Always evaluate integration quality quantitatively**
- **Compare multiple methods** when results are critical
- **Document your decision process**
- **Integration is a means to an end**, not the goal itself

**Most important:** Don't let integration remove real biological signal!
