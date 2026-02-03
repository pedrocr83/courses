# Annotation Confidence Framework

**A systematic approach to assigning and documenting cell type annotation confidence levels**

---

## Overview

Cell type annotation is **interpretation**, not fact. This framework helps you:
1. Systematically evaluate evidence for each annotation
2. Assign confidence levels based on multiple criteria
3. Document your reasoning
4. Communicate uncertainty appropriately

**Key Principle:** *It's better to report "likely CD4+ T cells" with honest uncertainty than to claim "CD4+ T cells" with false confidence.*

---

## Three-Tier Confidence System

### 🟢 **HIGH CONFIDENCE (Level 1)**

**Definition:** Strong, multiple lines of evidence support this annotation with minimal ambiguity.

**Criteria (Must meet at least 4):**
- [ ] ≥3 canonical markers clearly expressed
- [ ] Marker expression highly specific to this cluster (log2FC > 1.5)
- [ ] Markers absent in other clusters
- [ ] Matches published marker sets from same tissue
- [ ] Automated tools (e.g., SingleR) agree
- [ ] Cluster biologically plausible for tissue context
- [ ] Known negative markers are absent

**Example:**
```
Cluster 5 → "CD4+ T cells" [HIGH CONFIDENCE]

Evidence:
✓ CD3E+ CD4+ CD8- (canonical markers)
✓ IL7R+ LTB+ (additional T cell markers)
✓ CD8A/B negative (rules out CD8+ T cells)
✓ SingleR: T cells (score 0.89)
✓ Expected in PBMC sample
✓ log2FC (CD4) = 2.3 vs all other clusters
```

**Reporting:** Use definitive labels
- "CD4+ T cells"
- "Classical monocytes"
- "B cells"

---

### 🟡 **MEDIUM CONFIDENCE (Level 2)**

**Definition:** Reasonable evidence but with some ambiguity, missing markers, or borderline specificity.

**Criteria (Meets 2-3 of Level 1 criteria, OR):**
- [ ] Some canonical markers present but not all
- [ ] Markers present but not highly specific (log2FC 0.5-1.5)
- [ ] Automated tools partially agree or show moderate scores
- [ ] Represents a subtype that's harder to distinguish
- [ ] Literature support limited or conflicting

**Example:**
```
Cluster 8 → "Likely regulatory T cells (Tregs)" [MEDIUM CONFIDENCE]

Evidence:
✓ CD3E+ CD4+ (T cells confirmed)
✓ FOXP3+ (Treg marker, but low expression)
? IL2RA/CD25 low (expected to be high)
? No CTLA4 (another Treg marker)
✓ SingleR: T cells (score 0.72, subtype unclear)
⚠ Could be activated CD4+ T cells instead
```

**Reporting:** Use qualified labels
- "Likely regulatory T cells"
- "CD4+ T cells (possibly Tregs)"
- "Monocyte subset (CD16+?)"
- "Probable plasmacytoid DCs"

---

### 🔴 **LOW CONFIDENCE (Level 3)**

**Definition:** Weak evidence, high ambiguity, or novel population without clear markers.

**Criteria:**
- [ ] Few or no specific markers
- [ ] Markers conflict with expected biology
- [ ] Automated tools disagree or fail
- [ ] Could represent multiple cell types
- [ ] Potential doublets or low-quality cells
- [ ] Novel population not in literature

**Example:**
```
Cluster 12 → "Unknown myeloid-like cells" [LOW CONFIDENCE]

Evidence:
? LYZ+ (general myeloid marker)
? Mixed expression of monocyte and DC markers
✗ No clear canonical markers for any subtype
✗ SingleR: Low scores across multiple types (<0.5)
✗ High doublet score (0.35)
⚠ May be doublets or stressed cells
⚠ Consider sub-clustering or removing
```

**Reporting:** Use honest uncertainty
- "Unknown population"
- "Unidentified myeloid-like cells"
- "Cluster X (uncertain identity)"
- "Possible [cell type] (needs validation)"
- "Low-quality cells (consider filtering)"

---

## Evidence Scoring Matrix

Use this table to systematically score evidence for each annotation:

### Evidence Table Template

| Evidence Type | Score | Weight | Notes |
|--------------|-------|---------|-------|
| **Canonical markers (3+)** | 0-3 | ×2 | # markers clearly detected |
| **Marker specificity** | 0-3 | ×2 | log2FC relative to other clusters |
| **Negative markers** | 0-2 | ×1 | Expected absent markers are absent |
| **SingleR/automated** | 0-3 | ×1 | Score and consistency |
| **Literature support** | 0-2 | ×1 | Published marker sets match |
| **Biological plausibility** | 0-2 | ×1 | Makes sense in tissue context |
| **Sub-cluster coherence** | 0-2 | ×1 | Sub-clustering maintains identity |

**Total Score:** _____ / 30 possible

**Confidence Assignment:**
- **≥20 points** → HIGH CONFIDENCE (🟢)
- **12-19 points** → MEDIUM CONFIDENCE (🟡)
- **<12 points** → LOW CONFIDENCE (🔴)

---

## Detailed Evaluation Criteria

### 1. Canonical Markers (Most Important)

**What are canonical markers?**
- Well-established, widely cited markers for specific cell types
- Found in multiple publications and databases
- Protein-level validation exists (antibody staining, flow cytometry)

**How to score:**
- **3 points:** ≥3 canonical markers strongly expressed (mean expression > 1 in log-normalized)
- **2 points:** 2 canonical markers present
- **1 point:** 1 canonical marker present
- **0 points:** No canonical markers

**Examples:**

```
T cells:
- Essential: CD3E, CD3D
- Subtypes: CD4, CD8A, CD8B

B cells:
- Essential: CD79A, CD79B, MS4A1 (CD20)
- Plasma: MZB1, SDC1 (CD138)

Monocytes:
- Essential: LYZ, CD14, FCGR3A (CD16)
- Classical: CD14+ CD16-
- Non-classical: CD14low CD16+

NK cells:
- Essential: NKG7, GNLY, NCAM1 (CD56)
- Cytotoxic: GZMA, GZMB, PRF1
```

### 2. Marker Specificity

**How to assess:**

```python
# Calculate log2 fold change for marker in cluster vs all others
import scanpy as sc

# Find markers with minimum log2FC
sc.tl.rank_genes_groups(adata, 'cluster', method='wilcoxon', 
                         use_raw=False, 
                         pts=True)  # Include fraction detected

# Look at results
result = adata.uns['rank_genes_groups']
markers = sc.get.rank_genes_groups_df(adata, group='cluster_5')

# Check:
# - log2FC > 1.5 (strong specificity)
# - pct.1 > 50% (expressed in most cells of cluster)
# - pct.2 < 20% (rare in other clusters)
```

**Scoring:**
- **3 points:** log2FC > 2.0 for multiple markers, pct.1 > 70%, pct.2 < 10%
- **2 points:** log2FC > 1.5, pct.1 > 50%, pct.2 < 20%
- **1 point:** log2FC > 1.0, pct.1 > 30%, pct.2 < 30%
- **0 points:** log2FC < 1.0 or markers widely expressed

### 3. Negative Markers (Rules Out Other Types)

**Why important:** Confirming absence of markers for other cell types is as important as presence.

**Examples:**

```
Annotating CD4+ T cells?
Check that CD8A/CD8B are ABSENT (rules out CD8+ T cells)

Annotating CD8+ T cells?
Check that CD4 is ABSENT (rules out CD4+ T cells)

Annotating B cells?
Check that CD3E is ABSENT (rules out T cells)
Check that CD14 is ABSENT (rules out monocytes)

Annotating monocytes?
Check that CD3E is ABSENT (rules out T cells)
```

**Scoring:**
- **2 points:** All expected negative markers are clearly absent
- **1 point:** Most negative markers absent, some low expression
- **0 points:** Unexpected markers present (co-expression problem)

⚠️ **Warning:** If you see co-expression of lineage-incompatible markers (e.g., CD3E+ MS4A1+), suspect **doublets**!

### 4. Automated Annotation Agreement

**Tools to use:**
- **SingleR:** Reference-based annotation
- **CellTypist:** Marker-based annotation
- **Azimuth:** Reference mapping
- **scType:** Database-driven annotation

**How to interpret scores:**

```python
# SingleR example
import celltypist

# Load model and predict
model = celltypist.models.Model.load('Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model=model)

# Check scores
# High confidence: score > 0.8
# Medium: 0.5-0.8
# Low: < 0.5
```

**Scoring:**
- **3 points:** Multiple tools agree, scores > 0.8
- **2 points:** Tools partially agree, scores 0.6-0.8
- **1 point:** Tools agree weakly, scores 0.4-0.6
- **0 points:** Tools disagree or scores < 0.4

### 5. Literature Support

**Where to check:**
- **PanglaoDB:** https://panglaodb.se/ (cell type marker database)
- **CellMarker:** http://bio-bigdata.hrbmu.edu.cn/CellMarker/
- **PubMed:** Search "[cell type] markers [tissue]"
- **Original papers** describing similar datasets

**Scoring:**
- **2 points:** Multiple publications support exact marker combination
- **1 point:** General support but not exact match
- **0 points:** No literature support or contradicts literature

### 6. Biological Plausibility

**Ask yourself:**
- Is this cell type expected in this tissue?
- At this frequency?
- With these other cell types present?
- In this biological context?

**Examples:**

```
✅ Plausible:
- T cells in blood (common)
- Neurons in brain (expected)
- Hepatocytes in liver (abundant)

⚠️ Questionable:
- Neurons in blood (very rare, contamination?)
- Plasma cells at 30% in healthy PBMC (too high)
- Cardiomyocytes in lung (wrong tissue)

❌ Implausible:
- Photoreceptors in liver
- Red blood cells after filtering (should be removed)
- Cell types from different species (unless intentional)
```

**Scoring:**
- **2 points:** Highly plausible, expected in tissue
- **1 point:** Plausible but unexpected frequency
- **0 points:** Biologically implausible (investigate!)

---

## Annotation Workflow

### Step 1: Generate Marker Lists

```python
# Find top markers for each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Visualize top markers
sc.pl.rank_genes_groups(adata, n_genes=20)

# Get marker dataframe for specific cluster
markers = sc.get.rank_genes_groups_df(adata, group='0')
print(markers.head(20))
```

### Step 2: Check Canonical Markers

```python
# Create list of canonical markers to check
canonical_markers = {
    'T_cells': ['CD3E', 'CD3D', 'CD3G'],
    'B_cells': ['CD79A', 'CD79B', 'MS4A1'],
    'Monocytes': ['CD14', 'FCGR3A', 'LYZ'],
    'NK_cells': ['NKG7', 'GNLY', 'NCAM1'],
    'DCs': ['FCER1A', 'CST3', 'CLEC10A']
}

# Visualize on UMAP
for cell_type, markers in canonical_markers.items():
    sc.pl.umap(adata, color=markers, ncols=3, 
               title=[f"{m} ({cell_type})" for m in markers])
```

### Step 3: Run Automated Annotation

```python
# SingleR (R) or CellTypist (Python)
import celltypist

# Download model
celltypist.models.download_models()

# Predict
model = celltypist.models.Model.load('Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model=model)

# Add to adata
adata.obs['celltypist'] = predictions.predicted_labels

# Visualize
sc.pl.umap(adata, color=['leiden', 'celltypist'])
```

### Step 4: Manual Curation

For each cluster:

```python
# Create annotation template
cluster_id = 0

annotation_template = f"""
Cluster {cluster_id} Annotation
==============================

Top Markers:
{markers.head(10).to_string()}

Canonical Marker Check:
- CD3E: {check}
- CD79A: {check}
- [add more]

Automated Annotation:
- SingleR: {result}
- CellTypist: {result}

Proposed Label: _____________
Confidence: HIGH / MEDIUM / LOW
Reasoning: _____________
"""
```

### Step 5: Fill Evidence Matrix

| Cluster | Proposed Label | Evidence Score | Confidence | Notes |
|---------|----------------|----------------|------------|-------|
| 0 | CD4+ T cells | 24 | 🟢 HIGH | Clear CD3+CD4+CD8- |
| 1 | B cells | 22 | 🟢 HIGH | CD79A+ MS4A1+ CD3- |
| 2 | Classical monocytes | 20 | 🟢 HIGH | CD14+ CD16- LYZ+ |
| 3 | NK cells | 18 | 🟡 MEDIUM | NKG7+ but low GNLY |
| 4 | Likely CD8+ T cells | 15 | 🟡 MEDIUM | CD8A+ but low expression |
| 5 | Unknown myeloid | 8 | 🔴 LOW | Mixed markers, possible doublets |

---

## Special Cases

### Case 1: Doublet Clusters

**Symptoms:**
- Co-expression of lineage-incompatible markers (CD3E+ MS4A1+)
- Very high UMI/gene counts
- High doublet scores from Scrublet/DoubletFinder

**Action:**
```
Label: "Potential doublets (T/B cell)"
Confidence: N/A
Decision: REMOVE or mark for exclusion from downstream analysis
```

### Case 2: Low-Quality Clusters

**Symptoms:**
- Dominated by mitochondrial genes
- Very low marker specificity
- High ribosomal gene percentage
- Enriched for stress response genes

**Action:**
```
Label: "Low-quality cells"
Confidence: N/A
Decision: REMOVE or re-evaluate QC thresholds
```

### Case 3: Novel/Rare Populations

**Symptoms:**
- No clear canonical markers
- Not in reference databases
- Biologically plausible but undescribed

**Action:**
```
Label: "Unknown population (CD45+ LYZ+)"
Confidence: 🔴 LOW
Notes: Potential novel subtype
Recommendation: Further validation needed
- Isolate and re-sequence
- Protein validation
- Functional assays
```

### Case 4: Transitional States

**Symptoms:**
- Express markers from multiple states
- Fall between canonical types on UMAP
- Markers suggest differentiation/activation

**Action:**
```
Label: "Activated T cells (transitioning)"
Confidence: 🟡 MEDIUM
Notes: Intermediate activation state
- CD3E+ (T cell)
- CD69+ (activation)
- IL2RA+ (activation)
- Between naive and effector on UMAP
```

### Case 5: Cluster Should Be Split

**Symptoms:**
- Two sub-populations visible on UMAP
- Markers show heterogeneity within cluster
- Automated tools assign different labels

**Action:**
```
Recommendation: Sub-cluster!

# Subset cluster and re-cluster
cluster_3 = adata[adata.obs['leiden'] == '3']
sc.tl.leiden(cluster_3, resolution=0.5, key_added='sub_leiden')

# Annotate sub-clusters separately
```

### Case 6: Clusters Should Be Merged

**Symptoms:**
- No distinguishing markers between clusters
- Arbitrary split at current resolution
- Same cell type by all evidence

**Action:**
```
Clusters 2, 5, 8 → Merge to "CD8+ T cells"
Confidence: 🟢 HIGH (for merged identity)
Notes: Over-clustering artifact, no biological justification for split
```

---

## Documentation Template

### Annotation Report for [Dataset Name]

```markdown
# Cell Type Annotation Report

## Dataset Information
- Samples: ___
- Total cells (post-QC): ___
- Number of clusters: ___
- Tissue type: ___
- Sequencing platform: ___

## Annotation Summary

| Cluster ID | Cell Type | Confidence | % of Cells | Median Genes | Key Markers |
|------------|-----------|------------|------------|--------------|-------------|
| 0 | CD4+ T cells | 🟢 HIGH | 32% | 1,234 | CD3E, CD4, IL7R |
| 1 | B cells | 🟢 HIGH | 18% | 1,567 | CD79A, MS4A1, CD19 |
| ... | ... | ... | ... | ... | ... |

## High Confidence Annotations (Level 1)
[Detailed evidence for each]

## Medium Confidence Annotations (Level 2)
[Detailed evidence, alternative interpretations]

## Low Confidence Annotations (Level 3)
[Describe uncertainty, recommendations for follow-up]

## Decisions Made
- Clusters merged: [list]
- Clusters removed: [list] (reason: doublets/low-quality)
- Clusters requiring sub-clustering: [list]

## Validation Plan
- [ ] Check markers with protein data (if available)
- [ ] Compare to published datasets from same tissue
- [ ] Functional validation for [specific uncertain populations]

## Figures
1. UMAP colored by cluster + annotations
2. Dotplot of canonical markers
3. Heatmap of top markers per cluster
4. Confidence level distribution

## Caveats & Limitations
[List any concerns, ambiguities, or limitations]

## Conclusion
[Summary statement about annotation quality and readiness for downstream analysis]
```

---

## Best Practices

### DO:
✅ Start with broad categories, refine iteratively
✅ Check both positive AND negative markers
✅ Use multiple automated tools and compare
✅ Consult literature for your specific tissue
✅ Document evidence for EVERY annotation
✅ Report confidence levels honestly
✅ Revisit annotations after downstream analysis
✅ Get feedback from biologists familiar with the tissue

### DON'T:
❌ Copy labels from automated tools without validation
❌ Over-interpret low-confidence assignments
❌ Ignore co-expression of incompatible markers
❌ Claim high confidence without sufficient evidence
❌ Force every cluster into a canonical category
❌ Use annotation as final biological truth
❌ Skip documentation ("I'll remember")

---

## Confidence Level Communication Guide

### In Publications:

**HIGH CONFIDENCE:**
> "We identified CD4+ T cells (cluster 0) based on expression of CD3E, CD4, and IL7R, and absence of CD8A/B."

**MEDIUM CONFIDENCE:**
> "Cluster 3 likely represents regulatory T cells based on CD4 and FOXP3 expression, though canonical Treg markers CTLA4 and IL2RA showed moderate expression."

**LOW CONFIDENCE:**
> "Cluster 7 expressed myeloid markers LYZ and CST3 but lacked clear subtype-specific markers. This population requires further characterization."

### In Figures:

```
Use consistent visual indicators:
🟢 = HIGH confidence (solid labels)
🟡 = MEDIUM confidence (labels with "likely" or "?")
🔴 = LOW confidence (labels with "unknown" or "uncertain")
```

### In Supplementary Materials:

Include full evidence table with:
- All markers checked
- Automated tool results
- Confidence scores
- Alternative interpretations

---

## Quick Reference Checklist

For each cluster annotation:

- [ ] **Top 20 markers identified** and reviewed
- [ ] **≥3 canonical markers checked** (positive)
- [ ] **Negative markers checked** (rules out other types)
- [ ] **Automated annotation run** (SingleR/CellTypist)
- [ ] **Literature search performed** for tissue-specific markers
- [ ] **Biological plausibility assessed**
- [ ] **Evidence score calculated** (from matrix)
- [ ] **Confidence level assigned** (HIGH/MEDIUM/LOW)
- [ ] **Alternative interpretations considered**
- [ ] **Doublet potential evaluated**
- [ ] **Documentation completed**

---

## Summary

**Key Takeaways:**
1. **Annotation is interpretation**, not ground truth
2. **Confidence levels communicate uncertainty** appropriately
3. **Multiple lines of evidence** are essential
4. **Documentation** enables reproducibility and critique
5. **Honest uncertainty** is better than false confidence
6. **Validation** should guide future experiments

**Remember:** The goal is not perfect annotations, but **defensible annotations with clear evidence and appropriate confidence levels**.

---

**Next Steps:**
- Review `../labs/lab07_markers.ipynb` for marker identification
- See `../labs/lab08_manual_annotation.ipynb` for hands-on practice
- Check `../templates/annotation_report_template.md` for full report format
