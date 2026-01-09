# Course: Cell-Cell Communication Analysis

## From Biological Foundations to Computational Inference

---

## Course Overview

| | |
|---|---|
| **Duration** | 4 weeks (20-25 hours total) |
| **Format** | Self-paced with hands-on labs |
| **Level** | Intermediate to Advanced |
| **Prerequisites** | scRNA-seq processing, basic R/Python, clustering & annotation |

## Start Here (Do This First)

- **Step-by-step checklist**: `START_HERE.md`
- **Glossary**: `glossary.md`
- **Environment setup**: `setup.md`

### Learning Objectives

By the end of this course, you will be able to:
1. Understand the biological basis of cell-cell communication
2. Run and interpret CCC inference tools (CellChat, CellPhoneDB, NicheNet)
3. Critically compare results across methods
4. Build and visualize cell communication networks
5. Integrate spatial context into CCC analysis

---

## Course Structure

### Week 1: Biological Foundations & Data Preparation

#### Module 1: What is Cell-Cell Communication? (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 1.1 | Why cells talk: coordination, development, disease | Lecture |
| 1.2 | Signaling modes: juxtacrine, paracrine, autocrine, endocrine | Lecture |
| 1.3 | Key molecules: ligands, receptors, co-receptors | Lecture |
| 1.4 | Signaling cascades overview | Reading |
| **Lab 1** | Explore ligand-receptor databases | Hands-on |

**Key Concepts:** Signaling modalities, molecular players, biological context

---

#### Module 2: From Expression to Interaction (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 2.1 | The core assumption: co-expression implies communication | Lecture |
| 2.2 | Limitations: expression ≠ protein ≠ signaling | Lecture |
| 2.3 | Ligand-receptor databases: CellPhoneDB, CellChatDB, LIANA | Lecture |
| 2.4 | Statistical approaches to scoring interactions | Lecture |
| **Lab 2** | Prepare scRNA-seq data for CCC analysis | Hands-on |

**Key Concepts:** L-R pairing, database curation, statistical enrichment

---

#### Module 3: scRNA-seq Prerequisites for CCC (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 3.1 | Cell type annotation quality matters | Lecture |
| 3.2 | Handling dropout and low expression | Lecture |
| 3.3 | Choosing the right resolution | Discussion |
| **Lab 3** | Annotate cell types for CCC input | Hands-on |

**Key Concepts:** Annotation depth, expression thresholds, cluster granularity

---

### Week 2: Core CCC Inference Tools

#### Module 4: CellPhoneDB (2.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 4.1 | CellPhoneDB overview and database | Lecture |
| 4.2 | Statistical framework: permutation testing | Lecture |
| 4.3 | Running CellPhoneDB (Python) | Demo |
| 4.4 | Interpreting output: means, p-values, interactions | Lecture |
| **Lab 4** | Run CellPhoneDB on PBMC data | Hands-on |

**Key Concepts:** Permutation tests, interaction scores, dot plots

---

#### Module 5: CellChat (3 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 5.1 | CellChat philosophy: probability-based inference | Lecture |
| 5.2 | Communication probability and pathways | Lecture |
| 5.3 | Network centrality: senders, receivers, influencers | Lecture |
| 5.4 | Running CellChat (R) | Demo |
| 5.5 | Visualization: circle plots, hierarchy, heatmaps | Demo |
| **Lab 5A** | Run CellChat on PBMC data | Hands-on |
| **Lab 5B** | Compare conditions in CellChat | Hands-on |

**Key Concepts:** Communication probability, pathway aggregation, network roles

---

#### Module 6: NicheNet (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 6.1 | NicheNet's unique approach: ligand → target genes | Lecture |
| 6.2 | Prior knowledge integration: signaling + GRNs | Lecture |
| 6.3 | Ligand activity prediction | Lecture |
| 6.4 | Running NicheNet (R) | Demo |
| **Lab 6** | Predict ligand activities with NicheNet | Hands-on |

**Key Concepts:** Ligand-target links, prior networks, activity scores

---

### Week 3: Comparison, Integration & Spatial CCC

#### Module 7: LIANA - Consensus Framework (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 7.1 | Why consensus? The reproducibility problem | Lecture |
| 7.2 | LIANA's unified interface | Lecture |
| 7.3 | Aggregate scoring across methods | Lecture |
| **Lab 7** | Run LIANA and compare to individual tools | Hands-on |

**Key Concepts:** Method disagreement, consensus ranking, robustness

---

#### Module 8: Benchmarking & Critical Evaluation (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 8.1 | Why do tools give different answers? | Lecture |
| 8.2 | Benchmark studies overview | Reading |
| 8.3 | Choosing the right tool for your question | Discussion |
| **Lab 8** | Systematic comparison of 3 tools on same data | Hands-on |

**Key Concepts:** Method biases, gold standards, interpretation caveats

---

#### Module 9: Spatial Cell-Cell Communication (2.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 9.1 | Why spatial context matters | Lecture |
| 9.2 | Spatial technologies: Visium, MERFISH, Slide-seq | Lecture |
| 9.3 | Distance-weighted interaction scoring | Lecture |
| 9.4 | Spatial CCC tools overview | Lecture |
| **Lab 9** | Integrate spatial information into CCC | Hands-on |

**Key Concepts:** Proximity constraints, niche identification, gradients

---

### Week 4: Networks, ML & Applications

#### Module 10: CCC as a Graph Problem (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 10.1 | Graph formulation: nodes, edges, weights | Lecture |
| 10.2 | Network metrics: centrality, hubs, communities | Lecture |
| 10.3 | Visualization strategies | Demo |
| **Lab 10** | Build and analyze CCC networks | Hands-on |

**Key Concepts:** Graph construction, network analysis, visualization

---

#### Module 11: Machine Learning for CCC (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 11.1 | Embedding-based L-R scoring | Lecture |
| 11.2 | Graph neural networks for signaling | Lecture |
| 11.3 | Predictive models: perturbation outcomes | Lecture |
| **Lab 11** | Explore ML-enhanced CCC scoring | Hands-on |

**Key Concepts:** GNNs, embeddings, prediction vs inference

---

#### Module 12: Biological Interpretation & Applications (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 12.1 | From networks to hypotheses | Lecture |
| 12.2 | Disease contexts: cancer, inflammation, development | Case studies |
| 12.3 | Validation strategies | Lecture |
| 12.4 | Reporting CCC findings | Discussion |
| **Lab 12** | Interpret CCC in a disease context | Hands-on |

**Key Concepts:** Hypothesis generation, validation, publication standards

---

## Assessments

### Quizzes (after each module)
- 5-10 questions covering key concepts
- Immediate feedback

### Practical Assignments

| Assignment | Module | Description |
|------------|--------|-------------|
| **A1** | 2-3 | Prepare annotated data for CCC |
| **A2** | 4-5 | CellPhoneDB vs CellChat comparison |
| **A3** | 6-7 | NicheNet ligand activity + LIANA consensus |
| **A4** | 10 | Build and analyze a CCC network |

### Final Project
**Complete CCC analysis on a disease dataset**
- Choose cancer, inflammation, or development context
- Run 2+ methods
- Build interaction network
- Write biological interpretation (1-2 pages)

---

## Software Setup

### R Environment

```r
# Core CCC packages
install.packages("devtools")
devtools::install_github("sqjin/CellChat")
devtools::install_github("saeyslab/nichenetr")

# Dependencies
BiocManager::install(c("SingleCellExperiment", "ComplexHeatmap"))
install.packages(c("ggplot2", "igraph", "circlize"))
```

### Python Environment

```bash
pip install cellphonedb scanpy anndata pandas numpy
pip install liana  # Python version
```

### LIANA (R)

```r
devtools::install_github("saezlab/liana")
```

---

## Demo Datasets

| Dataset | Type | Source | Use |
|---------|------|--------|-----|
| PBMC 3k | scRNA-seq | 10x/Scanpy | Labs 1-8 |
| COVID PBMC | scRNA-seq | GEO | Disease context |
| Tumor microenvironment | scRNA-seq | Published study | Final project |

---

## Resources

### Primary References
- Armingol et al., 2021 - CCC inference review
- Jin et al., 2021 - CellChat paper
- Browaeys et al., 2020 - NicheNet paper
- Dimitrov et al., 2022 - LIANA paper

### Database Resources
| Resource | URL |
|----------|-----|
| CellPhoneDB | https://www.cellphonedb.org/ |
| CellChatDB | https://sqjin.github.io/CellChat/ |
| OmniPath | https://omnipathdb.org/ |

---

## Common Mistakes to Avoid

| Mistake | Consequence | Prevention |
|---------|-------------|------------|
| Poor cell annotation | Meaningless interactions | Validate annotations first |
| Ignoring low expression | Miss or hallucinate signals | Set appropriate thresholds |
| Single tool reliance | Biased results | Compare multiple methods |
| No biological validation | Overconfident conclusions | Check literature, design experiments |
| Ignoring spatial context | Miss proximity requirements | Consider tissue architecture |
| Over-interpreting networks | False discovery | Focus on high-confidence edges |

---

## What Bad CCC Analysis Looks Like (and How to Catch It)

### Symptoms (Red Flags)

- **Garbage-in**: CCC results change drastically when you slightly change clustering/annotation.
- **Hallucinated interactions**: top interactions involve genes barely expressed; interactions driven by tiny cell groups.
- **Tool disagreement**: methods disagree completely and you pick one without justification.
- **No validation story**: results are reported as facts without proposing validation or checking known biology.
- **Spatial mismatch**: strong “communication” predicted between cell types that are not co-localized in tissue.

### Common Failure Modes → What It Looks Like → Fix

- **Annotation quality too low / wrong granularity**  
  - **Looks like**: interactions between mislabeled cell types; “new” signals vanish after relabeling.  
  - **Fix**: revisit clustering/annotation; merge/split to biologically meaningful units; document confidence.

- **Expression thresholding errors**  
  - **Looks like**: interactions dominated by dropout noise (or missed due to too strict thresholds).  
  - **Fix**: tune min expression / percent-expressed cutoffs; check sensitivity; report thresholds.

- **Pseudoreplication (cells treated as replicates across samples)**  
  - **Looks like**: condition differences appear extremely significant, not reproducible across donors.  
  - **Fix**: use sample-aware comparisons; analyze per-sample or pseudobulk where possible.

- **Database mismatch / outdated curation**  
  - **Looks like**: interactions that conflict with known receptor-ligand biology; missing canonical pathways.  
  - **Fix**: compare databases (CellPhoneDB vs CellChatDB vs OmniPath); report database version.

- **Over-interpretation of network plots**  
  - **Looks like**: “hub” conclusions from visualization alone.  
  - **Fix**: quantify (centrality metrics), show robustness across methods, prioritize high-confidence edges.

### Minimum “Sanity Checks” Before You Report

- **Known biology**: at least a few canonical interactions/pathways are recovered (positive controls).
- **Robustness**: top interactions persist across reasonable parameters and/or across methods.
- **Localization**: if tissue context exists, interactions are plausible given proximity constraints.
- **Validation plan**: propose at least one orthogonal validation (protein, perturbation, imaging).

---

## Learning Path Progression

```
BEGINNER (Weeks 1-2)
│
├── Understand CCC biology
├── Run CellChat on demo data
├── Interpret basic outputs
└── Visualize interactions
        │
        ▼
INTERMEDIATE (Weeks 3-4)
│
├── Compare multiple tools
├── Use LIANA for consensus
├── Build CCC networks
└── Integrate spatial data
        │
        ▼
ADVANCED (Post-course)
│
├── ML-enhanced CCC
├── Multi-condition comparisons
├── Causal inference
└── Custom pipeline development
```

---

## Connection to Other Courses

| Prerequisite from | Why |
|-------------------|-----|
| scRNA-seq Processing | Need annotated, QC'd data |
| Differential Expression | DEGs define sender/receiver genes |

| Feeds into | Why |
|------------|-----|
| Trajectory Analysis | Signaling along differentiation |
| Drug Discovery | Target identification |
| Disease Modeling | Understand pathological communication |

---

## Suggested Schedule

| Day | Activity | Time |
|-----|----------|------|
| Mon/Wed/Fri | Lectures + Reading | 1 hr |
| Tue/Thu | Labs | 1.5 hrs |
| Weekend | Assignments | 2 hrs |

*Total: ~6 hours/week for 4 weeks*

---

## Certificate Criteria

To complete this course:
- [ ] Pass all module quizzes (≥70%)
- [ ] Submit all 4 practical assignments
- [ ] Complete final project with network and interpretation
- [ ] Document reproducible analysis

---

*This course bridges molecular biology and computational methods, preparing you to extract biological insight from cell communication networks.*

