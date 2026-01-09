# Assignment 1: Design Matrix Construction

**Module:** 3 - Experimental Design  
**Due:** End of Week 1  
**Points:** 100

---

## Overview

Construct appropriate design matrices for three experimental scenarios and justify your choices.

---

## Scenario A: Simple Two-Group Comparison (30 points)

**Study:** Researchers compare gene expression in tumor vs normal tissue from 6 patients.

| Sample | Tissue | Patient |
|--------|--------|---------|
| S1 | tumor | P1 |
| S2 | normal | P1 |
| S3 | tumor | P2 |
| S4 | normal | P2 |
| S5 | tumor | P3 |
| S6 | normal | P3 |
| ... | ... | ... |

**Tasks:**
1. Is this a paired design? Why?
2. Write the appropriate design formula for DESeq2
3. What is the reference level for the tissue factor?
4. How would you extract the tumor vs normal comparison?

---

## Scenario B: Batch Effect Correction (35 points)

**Study:** Drug treatment effect measured across two processing batches.

| Sample | Treatment | Batch |
|--------|-----------|-------|
| S1 | control | batch1 |
| S2 | control | batch1 |
| S3 | drug | batch1 |
| S4 | drug | batch1 |
| S5 | control | batch2 |
| S6 | control | batch2 |
| S7 | drug | batch2 |
| S8 | drug | batch2 |

**Tasks:**
1. Why is batch correction necessary here?
2. Write the design formula
3. Which factor should come first in the formula? Why?
4. Are there any confounding issues with this design?
5. Write the R code to set up the DESeqDataSet

---

## Scenario C: Factorial Design (35 points)

**Study:** Effect of drug treatment in wild-type (WT) and knockout (KO) mice.

| Sample | Treatment | Genotype |
|--------|-----------|----------|
| S1 | vehicle | WT |
| S2 | vehicle | WT |
| S3 | drug | WT |
| S4 | drug | WT |
| S5 | vehicle | KO |
| S6 | vehicle | KO |
| S7 | drug | KO |
| S8 | drug | KO |

**Tasks:**
1. What questions can this design answer?
2. Write the design formula for:
   a. Main effects only
   b. Main effects + interaction
3. Explain what the interaction term tests
4. Write the contrasts to extract:
   - Drug effect in WT
   - Drug effect in KO
   - Difference in drug effect between genotypes

---

## Deliverables

1. **Written answers** for all scenarios (PDF or Markdown)
2. **R code** demonstrating how to set up each DESeqDataSet

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Scenario A correct design | 15 |
| Scenario A justification | 15 |
| Scenario B correct design | 20 |
| Scenario B batch reasoning | 15 |
| Scenario C factorial design | 20 |
| Scenario C interaction interpretation | 15 |

---

## Hints

- Use `relevel()` to set reference levels
- Order matters: put nuisance variables before factors of interest
- Interaction terms use `:` or `*` in R formulas

