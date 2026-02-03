# Module 13: Cell-Cell Communication Validation & Experimental Design

**From computational predictions to validated biology**

---

## Overview

CCC inference tools predict thousands of interactions, but **prediction ≠ validation**. This module teaches you:
1. How to critically evaluate CCC predictions
2. Design validation experiments
3. Prioritize interactions for follow-up
4. Appropriately report CCC findings with caveats

**Duration:** 2-3 hours  
**Prerequisites:** Modules 4-7 (CCC inference tools)

---

## Learning Objectives

By the end of this module, you will be able to:
- [ ] Distinguish between high vs low-confidence CCC predictions
- [ ] Design validation experiments for predicted interactions
- [ ] Prioritize interactions based on multiple criteria
- [ ] Integrate CCC predictions with orthogonal data
- [ ] Report CCC results with appropriate caveats and validation plans

---

## The Validation Challenge

### Why CCC Validation is Hard

**CCC tools predict based on:**
- Co-expression of ligand (sender) and receptor (receiver)
- Database-curated ligand-receptor pairs
- Statistical enrichment or scoring

**But they DON'T directly measure:**
- ❌ Protein-level expression (using RNA as proxy)
- ❌ Protein secretion/presentation
- ❌ Receptor activation
- ❌ Downstream signaling
- ❌ Functional consequence
- ❌ Spatial proximity (in non-spatial datasets)

**Result:** High false positive rate without validation!

---

## Confidence Scoring Framework

### Three-Tier System for CCC Predictions

#### 🟢 HIGH CONFIDENCE
**Criteria (must meet ≥4):**
- [ ] Detected by ≥2 independent methods
- [ ] High expression of both L and R (not just present)
- [ ] Known biological context (literature support)
- [ ] Spatial proximity confirmed (if spatial data available)
- [ ] Consistent across samples/replicates
- [ ] Part of known signaling pathway
- [ ] Experimentally validated in literature (same tissue/context)

**Example:**
```
CXCL12 (mesenchymal cells) → CXCR4 (T cells) in lymph node

✓ Detected by CellPhoneDB, CellChat, LIANA
✓ Both highly expressed (>50% cells, high mean expression)
✓ Well-established chemokine signaling
✓ T cells and mesenchymal cells co-localize in tissue
✓ Consistent across all 6 samples
✓ CXCL12-CXCR4 axis extensively validated in literature
```

---

#### 🟡 MEDIUM CONFIDENCE
**Criteria:**
- Detected by 1 method OR weak signal by multiple
- Moderate expression levels
- Some literature support but not in exact context
- Spatial relationship plausible but not confirmed
- Variable across samples

**Example:**
```
CCL5 (NK cells) → CCR5 (T cells)

✓ Detected by CellPhoneDB only
? Expression moderate (30% of cells)
✓ CCL5-CCR5 known interaction
? Not specifically validated in this tissue type
? Present in 4/6 samples (variable)
```

---

#### 🔴 LOW CONFIDENCE
**Criteria:**
- Weak signal or single method detection
- Low expression of L or R
- No literature support
- Biologically implausible
- Batch-specific or sample-specific

**Example:**
```
GENE_X (Cell Type A) → GENE_Y (Cell Type B)

✗ Only detected by one method, low score
✗ Low expression (5% of cells)
✗ No known interaction in literature
✗ Cell types not co-localized in tissue
✗ Only in 1 sample (likely artifact)
```

---

## Validation Strategy Pyramid

### Level 1: Computational Validation (Quick, Low Cost)

**Can be done with existing data**

#### 1.1 Method Consensus
```python
# Compare results across tools
from liana import Method

# Run multiple methods via LIANA
import liana as li

li.mt.rank_aggregate(
    adata,
    groupby='celltype',
    use_raw=False,
    verbose=True,
    methods=['cellphonedb', 'cellchat', 'natmi', 'connectome']
)

# Get consensus interactions
consensus = adata.uns['liana_res']
high_conf = consensus[consensus['magnitude_rank'] < 0.2]  # Top 20%

print(f"High confidence interactions: {len(high_conf)}")
```

---

#### 1.2 Expression Level Check
```python
import scanpy as sc

# For a candidate interaction: CXCL12 (Fibroblasts) → CXCR4 (T cells)
ligand = 'CXCL12'
receptor = 'CXCR4'
sender = 'Fibroblast'
receiver = 'T_cell'

# Check expression levels
sender_expr = adata[adata.obs['celltype'] == sender, ligand].X.toarray()
receiver_expr = adata[adata.obs['celltype'] == receiver, receptor].X.toarray()

print(f"{ligand} in {sender}:")
print(f"  Mean: {sender_expr.mean():.3f}")
print(f"  % expressing: {(sender_expr > 0).mean() * 100:.1f}%")

print(f"{receptor} in {receiver}:")
print(f"  Mean: {receiver_expr.mean():.3f}")
print(f"  % expressing: {(receiver_expr > 0).mean() * 100:.1f}%")

# Thresholds for high confidence:
# - Mean expression > 0.5 (log-normalized)
# - % expressing > 20%
```

---

#### 1.3 Literature Mining
```python
# Check if interaction is known
from Bio import Entrez

def check_pubmed(ligand, receptor, tissue=""):
    Entrez.email = "your_email@example.com"
    
    query = f"{ligand} AND {receptor} AND (signaling OR communication)"
    if tissue:
        query += f" AND {tissue}"
    
    handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
    record = Entrez.read(handle)
    handle.close()
    
    return int(record["Count"])

# Check your interaction
pubmed_count = check_pubmed("CXCL12", "CXCR4", "lymphoid")
print(f"PubMed papers: {pubmed_count}")

# High confidence: >10 papers
# Medium: 1-10 papers
# Low: 0 papers
```

---

#### 1.4 Pathway Coherence
```python
# Check if L-R pair is part of known pathway
import gseapy as gp

# Get downstream targets of receptor
# (from NicheNet, KEGG, or other sources)
targets = ['GENE1', 'GENE2', 'GENE3', ...]  # Known targets of CXCR4

# Check if targets are expressed in receiver cells
receiver_cells = adata[adata.obs['celltype'] == receiver]

target_expr = []
for gene in targets:
    if gene in receiver_cells.var_names:
        expr = receiver_cells[:, gene].X.toarray().mean()
        target_expr.append(expr)

print(f"Mean target expression: {np.mean(target_expr):.3f}")

# High confidence: Targets are expressed (pathway active)
# Low confidence: Targets not expressed (no downstream effect)
```

---

### Level 2: Orthogonal Data Integration (Moderate Cost)

**Requires additional data types**

#### 2.1 Spatial Transcriptomics
```python
# If you have spatial data, check co-localization
import squidpy as sq

# Load spatial data
adata_spatial = sc.read_visium('spatial_data/')

# Calculate spatial neighbors
sq.gr.spatial_neighbors(adata_spatial)

# Check if sender and receiver are spatially proximal
sender_spots = adata_spatial.obs['celltype'] == sender
receiver_spots = adata_spatial.obs['celltype'] == receiver

# Calculate spatial proximity score
# (Implementation depends on your spatial data format)

# High confidence: Cell types are spatial neighbors
# Low confidence: Cell types are distant
```

---

#### 2.2 Protein-Level Validation
```python
# If you have CITE-seq (protein) data
# Check if ligand/receptor expressed at protein level

if 'protein' in adata.obsm.keys():
    protein_data = adata.obsm['protein']
    
    # Check for matching antibody
    if 'CXCR4' in protein_data.columns:
        rna_expr = adata[:, 'CXCR4'].X.toarray()
        protein_expr = protein_data['CXCR4'].values
        
        from scipy.stats import spearmanr
        corr, pval = spearmanr(rna_expr, protein_expr)
        
        print(f"RNA-Protein correlation: r={corr:.3f}, p={pval:.2e}")
        
        # High confidence: Strong RNA-protein correlation (r > 0.5)
        # Low confidence: Poor correlation (protein absent despite RNA)
```

---

#### 2.3 Perturbation Data
```python
# If you have knockout or inhibitor data
# Check if blocking ligand affects receiver

# Compare receiver cells: control vs ligand-blocked
control = adata[adata.obs['condition'] == 'control']
blocked = adata[adata.obs['condition'] == 'CXCL12_blocked']

# Get receiver cell type
control_receiver = control[control.obs['celltype'] == receiver]
blocked_receiver = blocked[blocked.obs['celltype'] == receiver]

# Check downstream target expression
target = 'TARGET_GENE'  # Known downstream gene
control_expr = control_receiver[:, target].X.toarray().mean()
blocked_expr = blocked_receiver[:, target].X.toarray().mean()

print(f"Control: {control_expr:.3f}, Blocked: {blocked_expr:.3f}")

# High confidence: Blocking reduces target expression
# Low confidence: No change (interaction doesn't affect targets)
```

---

### Level 3: Experimental Validation (Gold Standard)

**Requires new experiments**

#### 3.1 Co-Culture Assays

**Design:**
1. Isolate sender and receiver cell populations
2. Culture separately (control)
3. Co-culture sender + receiver
4. Measure receiver response

**Readouts:**
- Receptor phosphorylation (Western blot)
- Downstream gene expression (qPCR, RNA-seq)
- Phenotypic changes (proliferation, migration, differentiation)

**Expected:**
- Control: Low baseline response
- Co-culture: Strong response
- Co-culture + blocking antibody: Response reduced

**Example Protocol:**
```
Sender cells (Fibroblasts) --[secrete CXCL12]--> Receiver cells (T cells)

Setup:
1. T cells alone (control)
2. T cells + Fibroblasts (co-culture)
3. T cells + Fibroblasts + anti-CXCL12 (blocking)

Measure: T cell migration (transwell assay) at 24h

Expected:
- Control: Minimal migration
- Co-culture: Increased migration (5-fold)
- Blocking: Migration reduced to control levels
```

---

#### 3.2 Ligand Recombinant Protein

**Design:**
1. Treat receiver cells with recombinant ligand
2. Measure receptor activation and downstream effects

**Advantages:**
- Direct test of ligand-receptor interaction
- Dose-response curves
- No need for sender cells

**Example:**
```
Treat T cells with recombinant CXCL12:
- 0 ng/mL (control)
- 10 ng/mL
- 100 ng/mL
- 1000 ng/mL

Measure:
- CXCR4 phosphorylation (30 min)
- ERK activation (Western blot)
- Target gene expression (2h, 6h, 24h)
- Cell migration (24h)

Expected:
- Dose-dependent increase in all readouts
- EC50 should match known values (~10-100 ng/mL for CXCL12)
```

---

#### 3.3 Receptor Blockade

**Design:**
1. Block receptor with antibody or small molecule
2. Observe effect on receiver cells in presence of sender

**Example:**
```
Experimental groups:
1. T cells + Fibroblasts (control)
2. T cells + Fibroblasts + anti-CXCR4 (blocking Ab)
3. T cells + Fibroblasts + isotype control (negative control)
4. T cells + Fibroblasts + AMD3100 (CXCR4 antagonist)

Measure:
- T cell activation markers (CD69, CD25)
- Proliferation (CFSE dilution)
- Cytokine secretion (ELISA)

Expected:
- Control: Activation
- Anti-CXCR4 / AMD3100: Reduced activation
- Isotype control: No change (same as control)
```

---

#### 3.4 Genetic Perturbation

**Design:**
1. Knockout or knockdown ligand in sender cells
2. Knockout or knockdown receptor in receiver cells
3. Compare to wild-type

**Example:**
```
CRISPR knockout:
- Sender: CXCL12 KO fibroblasts
- Receiver: CXCR4 KO T cells

Setup:
1. WT fibroblasts + WT T cells (control)
2. CXCL12 KO fibroblasts + WT T cells
3. WT fibroblasts + CXCR4 KO T cells
4. CXCL12 KO + CXCR4 KO (double negative control)

Measure: T cell migration, activation, function

Expected:
- Groups 2, 3, 4: Reduced response vs control
- Confirms both ligand and receptor are required
```

---

#### 3.5 Spatial Validation

**Design:**
1. Immunofluorescence or immunohistochemistry
2. Co-stain for cell type markers + ligand + receptor
3. Confirm co-localization in tissue

**Example:**
```
Staining panel:
- Cell type 1 marker: CD3 (T cells)
- Cell type 2 marker: Vimentin (Fibroblasts)
- Ligand: CXCL12
- Receptor: CXCR4

Analysis:
- Segment cells by type
- Measure distance between CXCL12+ fibroblasts and CXCR4+ T cells
- Check for co-localization (<50 μm distance)

Expected:
- T cells near CXCL12+ fibroblasts express higher CXCR4
- Distance-dependent relationship
```

---

## Prioritization Framework

**You can't validate everything! How to choose interactions for follow-up:**

### Priority Score Calculation

```python
import pandas as pd

# Example scoring system
def calculate_priority_score(interaction_row):
    score = 0
    
    # Computational evidence (0-30 points)
    if interaction_row['n_methods'] >= 3:
        score += 15
    elif interaction_row['n_methods'] == 2:
        score += 10
    else:
        score += 5
    
    score += min(interaction_row['mean_expression'] * 10, 15)  # Max 15
    
    # Biological relevance (0-40 points)
    if interaction_row['pubmed_papers'] > 10:
        score += 20
    elif interaction_row['pubmed_papers'] > 0:
        score += 10
    
    if interaction_row['known_pathway']:
        score += 10
    
    if interaction_row['disease_relevant']:
        score += 10
    
    # Feasibility (0-30 points)
    if interaction_row['protein_available']:  # Recombinant protein
        score += 10
    
    if interaction_row['antibody_available']:  # Blocking Ab
        score += 10
    
    if interaction_row['cell_type_isolatable']:  # Can isolate both cell types
        score += 10
    
    return score

# Apply to all predictions
interactions['priority_score'] = interactions.apply(calculate_priority_score, axis=1)

# Sort and prioritize
top_interactions = interactions.nlargest(10, 'priority_score')
```

---

### Decision Matrix

| Criterion | Weight | How to Score |
|-----------|--------|-------------|
| **Method consensus** | 15% | # of methods detecting |
| **Expression level** | 15% | Mean expr of L and R |
| **Literature support** | 20% | PubMed papers |
| **Known pathway** | 10% | In KEGG/Reactome? |
| **Disease relevance** | 10% | Relevant to research question |
| **Feasibility** | 30% | Tools available? Cells isolatable? |

**Total: 100%**

High priority: Score ≥ 70  
Medium priority: Score 50-69  
Low priority: Score < 50

---

## Validation Checklist

### Before Claiming an Interaction is "Validated"

- [ ] **Computational confidence**: Detected by ≥2 methods, high expression
- [ ] **Literature support**: ≥1 paper describing interaction in relevant context
- [ ] **Spatial proximity**: Co-localization confirmed (if spatial data available)
- [ ] **Expression validation**: Protein-level or independent dataset confirmation
- [ ] **Functional validation**: At least ONE of:
  - [ ] Co-culture + blocking experiment
  - [ ] Recombinant ligand treatment
  - [ ] Receptor blockade
  - [ ] Genetic perturbation
- [ ] **Downstream effects**: Measured pathway activation or phenotypic change
- [ ] **Reproducibility**: Validated in ≥2 independent experiments

**Minimum for publication:**
- ✅ Computational confidence + Literature support
- ✅ At least 1 functional validation
- ✅ Reproducibility across experiments

---

## Reporting CCC Results

### Publication Guidelines

#### Figure Legends

**For computational predictions (no validation):**
> "Predicted cell-cell interactions between [cell types] inferred using CellPhoneDB and CellChat. Dot size represents significance; color represents mean expression. Top 20 interactions shown based on consensus ranking. These are **computational predictions** and require experimental validation."

**With validation:**
> "Cell-cell communication network inferred from scRNA-seq data. Selected interactions (marked with asterisks) were validated by [method]. CXCL12-CXCR4 interaction between fibroblasts and T cells was confirmed by co-culture assays and recombinant protein stimulation (Fig. Xb)."

---

#### Methods Section

Include:
1. **CCC tools used** (with versions, parameters)
2. **Filtering criteria** (expression thresholds, p-value cutoffs)
3. **Prioritization method** (how top interactions selected)
4. **Validation methods** (if applicable)

**Example:**
> "Cell-cell communication was inferred using CellPhoneDB (v4.0) and CellChat (v1.6). Only interactions with mean expression > 0.5 in both sender and receiver and FDR < 0.05 were retained. Interactions were prioritized based on method consensus (detected by both tools) and literature support (PubMed search). The top 3 prioritized interactions were experimentally validated using co-culture assays and receptor blockade experiments (see Supplementary Methods)."

---

#### Results Section

**Be honest about confidence:**

🟢 **High confidence (validated):**
> "We identified CXCL12-CXCR4 signaling from fibroblasts to T cells. This interaction was robustly detected (CellPhoneDB p=0.001, CellChat probability=0.85) and validated experimentally. Co-culture of fibroblasts with T cells induced CXCR4 phosphorylation and downstream ERK activation (Fig. 3c), which was abolished by anti-CXCL12 blocking antibody (Fig. 3d)."

🟡 **Medium confidence (predicted, partially validated):**
> "Our analysis suggests CCL5-CCR5 interaction between NK cells and T cells. This was predicted by CellPhoneDB (p=0.03) and supported by literature in similar contexts. However, expression levels were moderate and experimental validation is needed to confirm functional significance."

🔴 **Low confidence (prediction only):**
> "Several additional interactions were predicted but showed low consensus across methods or limited expression. These remain candidate interactions requiring further investigation (Supplementary Table X)."

---

### Supplementary Materials

**Required:**
1. **Full interaction tables** (all predictions, not just top hits)
2. **Method comparison** (overlap, discrepancies)
3. **Expression heatmaps** for top interactions
4. **Validation data** (full experimental results)
5. **Prioritization criteria** and scores

---

## Common Pitfalls & How to Avoid

### Pitfall 1: Over-Interpreting Predictions

**Problem:**
> "Our analysis revealed that fibroblasts communicate with T cells via CXCL12-CXCR4 signaling."

**Issue:** "Revealed" implies validation. It's a prediction!

**Better:**
> "Our analysis **predicts** CXCL12-CXCR4 signaling between fibroblasts and T cells, supported by high expression of both ligand and receptor and detection by multiple methods."

---

### Pitfall 2: Ignoring Low Expression

**Problem:** Reporting interaction with ligand expressed in 2% of sender cells

**Solution:**
```python
# Filter by expression before reporting
min_pct = 10  # At least 10% of cells
min_expr = 0.5  # Mean expression > 0.5

interactions_filtered = interactions[
    (interactions['sender_pct'] > min_pct) &
    (interactions['receiver_pct'] > min_pct) &
    (interactions['mean_expr'] > min_expr)
]
```

---

### Pitfall 3: Cherry-Picking Methods

**Problem:** Running 5 methods, only reporting the one that found your favorite interaction

**Solution:**
- Report all methods used
- Show method agreement/disagreement
- Focus on consensus interactions
- Disclose if specific interaction only found by one method

---

### Pitfall 4: No Validation Plan

**Problem:** Presenting predictions without proposing how to validate

**Solution:**
Always include:
> "Validation plan: We propose to validate the top 3 predicted interactions using: (1) co-culture assays with blocking antibodies, (2) recombinant protein stimulation, and (3) spatial validation using multiplexed immunofluorescence. Preliminary data suggests [if you have any]. Full validation is planned as follow-up work."

---

## Validation Workflow Template

### Step 1: Generate Predictions
```python
# Run CCC analysis
# Filter to high-confidence
# Rank by priority score
```

### Step 2: Computational Validation
```python
# Check method consensus
# Verify expression levels
# Literature mining
# Pathway coherence
```

### Step 3: Prioritize Top 5-10
```python
# Calculate priority scores
# Consider biological relevance
# Assess feasibility
```

### Step 4: Integrate Orthogonal Data (if available)
```python
# Spatial data
# Protein data (CITE-seq)
# Perturbation data
```

### Step 5: Design Experiments
```python
# Choose 2-3 top interactions
# Plan validation strategy
# Order reagents
# Write protocol
```

### Step 6: Validate
```python
# Run experiments
# Collect data
# Analyze results
# Confirm or refute prediction
```

### Step 7: Report
```python
# Validated: High confidence claims
# Predicted only: Report with caveats
# Failed validation: Report negative results (important!)
```

---

## Resources

### Databases
- **CellPhoneDB:** Curated L-R pairs with evidence
- **CellChatDB:** Literature-supported interactions
- **KEGG:** Pathway context
- **Reactome:** Signaling pathways

### Reagents
- **Recombinant proteins:** R&D Systems, PeproTech
- **Blocking antibodies:** BioLegend, BD Biosciences
- **Small molecule inhibitors:** Tocris, Selleckchem

### Validation Examples
- Check supplementary materials of recent CCC papers
- scRNA-seq + validation: Nature, Cell, Immunity journals

---

## Summary

**Key Takeaways:**
1. **CCC predictions are hypotheses**, not facts
2. **Validation is essential** for confident biological claims
3. **Prioritize ruthlessly** - can't validate everything
4. **Multiple lines of evidence** increase confidence
5. **Report honestly** - distinguish predicted from validated
6. **Negative results matter** - not all predictions are real

**Golden Rule:**
> "A validated interaction in one system is worth 100 predicted interactions without validation."

---

**Next Steps:**
- Complete prioritization for your CCC analysis
- Design validation experiments for top 3 interactions
- See examples in `../assignments/assignment4_network_analysis.md`
- Review validation section in final project guidelines

**Remember:** The goal is biological insight, not a long list of predictions!
