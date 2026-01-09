# Design Matrix Template

## Experiment Information

| Field | Value |
|-------|-------|
| Experiment Name | |
| Date | |
| Analyst | |

---

## Sample Metadata

| sample_id | condition | batch | donor | sex | age | other |
|-----------|-----------|-------|-------|-----|-----|-------|
| | | | | | | |
| | | | | | | |
| | | | | | | |

---

## Design Considerations

### Primary Comparison
- **Factor of interest:** 
- **Reference level:** 
- **Test level(s):** 

### Covariates to Include
- [ ] Batch
- [ ] Donor/Subject (for paired designs)
- [ ] Sex
- [ ] Age
- [ ] Other: ___________

### Design Type
- [ ] Simple two-group comparison
- [ ] Multi-group comparison
- [ ] Paired/blocked design
- [ ] Factorial design (multiple factors)
- [ ] Time series

---

## Design Formula

### R Formula Syntax

```r
# Simple two-group
design <- ~ condition

# With batch effect
design <- ~ batch + condition

# Paired design (donor as blocking factor)
design <- ~ donor + condition

# Factorial design
design <- ~ factor_A * factor_B

# Time series with subject
design <- ~ subject + timepoint
```

### Selected Design

```r
design <- ~ 
```

**Justification:**

[Explain why this design was chosen]

---

## Contrast Definition

### DESeq2 Contrasts

```r
# Extract results for specific comparison
results(dds, contrast = c("condition", "treatment", "control"))

# Or using name
results(dds, name = "condition_treatment_vs_control")
```

### edgeR Contrasts

```r
# Define contrast
contrast <- makeContrasts(treatment - control, levels = design)
```

---

## Checklist

- [ ] All samples have metadata
- [ ] Reference levels are set correctly
- [ ] Batch effects are accounted for
- [ ] Sample sizes are balanced (or noted)
- [ ] Design formula is valid
- [ ] Contrasts are correctly specified

---

## Common Design Patterns

### 1. Simple Comparison

```
Condition A (n=3) vs Condition B (n=3)
Design: ~ condition
```

### 2. With Batch Correction

```
Condition A vs B, processed in 2 batches
Design: ~ batch + condition
```

### 3. Paired Samples

```
Before vs After treatment, same subjects
Design: ~ subject + treatment
```

### 4. Two Factors

```
Treatment (drug/control) × Genotype (WT/KO)
Design: ~ treatment * genotype
```

---

## Notes

[Additional notes about the experimental design]

