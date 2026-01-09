# Cell Type Annotation Report

## Dataset Information

| Property | Value |
|----------|-------|
| Dataset | |
| Total cells | |
| Total genes | |
| After QC cells | |
| Date | |

---

## Preprocessing Summary

| Step | Parameters | Result |
|------|------------|--------|
| QC filtering | min_genes=, max_genes=, max_mt%= | cells removed |
| Normalization | | |
| HVG selection | n_hvgs= | genes selected |
| PCA | n_pcs= | |

---

## Clustering Summary

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| n_neighbors | | |
| Resolution | | |
| Algorithm | Leiden/Louvain | |
| Clusters | | |

---

## Cell Type Annotations

| Cluster | Cells | % | Cell Type | Confidence | Key Markers | Evidence |
|---------|-------|---|-----------|------------|-------------|----------|
| 0 | | | | High/Med/Low | | |
| 1 | | | | | | |
| 2 | | | | | | |
| 3 | | | | | | |
| 4 | | | | | | |

---

## Marker Evidence

### Cluster X: [Cell Type Name]

**Top markers:**
1. GENE1 - log2FC=, pval=
2. GENE2 - log2FC=, pval=
3. GENE3 - log2FC=, pval=

**Known markers detected:**
- MARKER1: expressed (✓)
- MARKER2: expressed (✓)

**Literature support:**
- [Citation]

---

## Automated Annotation Comparison

| Cluster | Manual | SingleR | CellTypist | Final |
|---------|--------|---------|------------|-------|
| 0 | | | | |
| 1 | | | | |

Agreement rate: __%

---

## Ambiguous Clusters

### Cluster X

**Issue:** 

**Markers detected:**

**Resolution:**

---

## Quality Concerns

| Issue | Clusters Affected | Action |
|-------|------------------|--------|
| Possible doublets | | |
| Low confidence | | |
| Over-clustering | | |

---

## Final UMAP

[Insert figure]

---

## Recommendations

1. 
2. 
3. 

---

## Files Generated

- `adata_annotated.h5ad`
- `markers.csv`
- `annotations.csv`
- Figures in `/figures/`

