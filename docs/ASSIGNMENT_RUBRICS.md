# Assignment Rubrics

**Standardized grading criteria for all practical assignments across the curriculum**

---

## General Rubric Structure

All assignments are graded on 5 dimensions:

| Dimension | Weight | Description |
|-----------|--------|-------------|
| **Code Quality** | 20% | Readability, organization, efficiency |
| **Documentation** | 20% | Comments, explanations, reproducibility |
| **Correctness** | 30% | Accurate execution, appropriate methods |
| **Interpretation** | 20% | Biological insight, critical thinking |
| **Visualization** | 10% | Figure quality, clarity, completeness |

**Total:** 100 points

---

## Detailed Rubric by Dimension

### 1. Code Quality (20 points)

| Score | Criteria |
|-------|----------|
| **18-20 (Excellent)** | Clean, well-organized code with clear variable names; functions used appropriately; efficient implementations; follows style guidelines (PEP8/tidyverse); no redundant code |
| **15-17 (Good)** | Generally clean code with minor organizational issues; mostly clear variable names; some inefficiencies but functional; mostly follows style guidelines |
| **12-14 (Satisfactory)** | Code works but difficult to follow; inconsistent naming; some redundancy; limited use of functions; style guidelines partially followed |
| **8-11 (Needs Improvement)** | Poorly organized code; unclear variable names; significant redundancy; hard to understand logic flow; style guidelines ignored |
| **0-7 (Unsatisfactory)** | Code is unreadable, disorganized, or contains major structural problems |

**Key Elements to Check:**
- [ ] Meaningful variable names (not `x`, `y`, `temp`)
- [ ] Consistent formatting and indentation
- [ ] Appropriate use of functions/modularity
- [ ] No hard-coded values that should be parameters
- [ ] Efficient implementations (no unnecessary loops)
- [ ] Code runs without errors

---

### 2. Documentation (20 points)

| Score | Criteria |
|-------|----------|
| **18-20 (Excellent)** | Comprehensive comments explaining logic; clear README or notebook markdown; all parameters documented; software versions recorded; data sources cited; fully reproducible from provided files |
| **15-17 (Good)** | Good comments on key sections; README present; most parameters documented; versions mostly recorded; largely reproducible |
| **12-14 (Satisfactory)** | Basic comments; minimal documentation; some parameters explained; versions partially recorded; reproducibility requires some guesswork |
| **8-11 (Needs Improvement)** | Sparse comments; insufficient documentation; parameters not explained; versions missing; difficult to reproduce |
| **0-7 (Unsatisfactory)** | No meaningful comments; no documentation; not reproducible |

**Key Elements to Check:**
- [ ] Comments explain **why**, not just what
- [ ] Markdown cells (notebooks) or README explain workflow
- [ ] All packages/libraries listed with versions
- [ ] Data sources and download instructions provided
- [ ] Parameters and thresholds justified
- [ ] Reproducible: Can someone else run your code?

---

### 3. Correctness (30 points)

| Score | Criteria |
|-------|----------|
| **27-30 (Excellent)** | All analyses executed correctly; appropriate methods chosen; proper statistical tests applied; results accurate; handles edge cases; validates outputs |
| **23-26 (Good)** | Analyses mostly correct; appropriate methods; minor errors that don't affect conclusions; results largely accurate |
| **18-22 (Satisfactory)** | Basic analyses correct; some inappropriate method choices; several errors but main results sound; limited validation |
| **12-17 (Needs Improvement)** | Significant errors in execution; inappropriate methods; results questionable; little to no validation |
| **0-11 (Unsatisfactory)** | Major errors; incorrect methods; results invalid or analysis incomplete |

**Key Elements to Check:**
- [ ] Correct implementation of required analyses
- [ ] Appropriate choice of methods for data type
- [ ] Proper use of statistical tests (when applicable)
- [ ] Results match expectations or explained if not
- [ ] Quality control applied appropriately
- [ ] Output files in correct format

---

### 4. Interpretation (20 points)

| Score | Criteria |
|-------|----------|
| **18-20 (Excellent)** | Insightful biological interpretation; connects results to broader context; identifies limitations; proposes follow-up experiments; demonstrates critical thinking; caveats clearly stated |
| **15-17 (Good)** | Solid interpretation; biological context provided; some limitations noted; reasonable conclusions; mostly appropriate caveats |
| **12-14 (Satisfactory)** | Basic interpretation present; limited biological context; few limitations noted; conclusions stated but not deeply justified |
| **8-11 (Needs Improvement)** | Superficial interpretation; little biological context; no limitations discussed; conclusions questionable or unsupported |
| **0-7 (Unsatisfactory)** | No interpretation or interpretation contradicts results; biological context absent |

**Key Elements to Check:**
- [ ] Results explained in biological terms (not just statistics)
- [ ] Connections to known biology/literature
- [ ] Limitations and caveats acknowledged
- [ ] Alternative explanations considered
- [ ] Conclusions supported by results
- [ ] Appropriate confidence level in claims

---

### 5. Visualization (10 points)

| Score | Criteria |
|-------|----------|
| **9-10 (Excellent)** | Clear, publication-quality figures; axes labeled with units; legends present; color schemes appropriate; multiple figure types used effectively; figures directly support conclusions |
| **7-8 (Good)** | Clear figures; mostly labeled; legends present; appropriate color schemes; figures support analysis |
| **5-6 (Satisfactory)** | Basic figures present; some labeling; basic legends; default colors; figures show results but not optimized |
| **3-4 (Needs Improvement)** | Unclear figures; poor labeling; legends missing or unclear; inappropriate colors; figures don't effectively communicate |
| **0-2 (Unsatisfactory)** | Figures missing, unreadable, or misleading |

**Key Elements to Check:**
- [ ] All axes labeled with variable names and units
- [ ] Titles informative
- [ ] Legends present and clear
- [ ] Color schemes accessible (colorblind-friendly when possible)
- [ ] Figure size appropriate (readable)
- [ ] Multiple complementary visualization types
- [ ] Figures referenced in text

---

## Course-Specific Rubrics

### Course 1: Single-Cell Processing

#### Assignment 1: FASTQ Inspection Report

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **FASTQ Structure** | 15 | Correctly identifies read structure (barcode, UMI, cDNA); shows example reads; describes quality scores |
| **Read Quality** | 15 | Analyzes quality score distributions; identifies any quality issues; proposes filtering if needed |
| **Chemistry Identification** | 10 | Correctly identifies 10x chemistry version or protocol; justifies conclusion |
| **Barcode Analysis** | 15 | Examines barcode diversity; checks against whitelist; identifies potential issues |
| **Code Quality** | 15 | (Standard rubric) |
| **Documentation** | 15 | (Standard rubric) |
| **Visualization** | 10 | Quality score plots; barcode rank plots |
| **Interpretation** | 5 | Assesses overall data quality; recommendations for proceeding |

---

#### Assignment 2: Pipeline Comparison

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Cell Ranger Results** | 15 | Correctly runs Cell Ranger; reports key metrics; saves outputs |
| **kb Results** | 15 | Correctly runs kallisto\|bustools; reports comparable metrics |
| **Comparison Table** | 20 | Side-by-side comparison of: cell count, gene count, mapping rate, runtime, memory usage |
| **Agreement Analysis** | 15 | Compares count matrices; calculates correlation; identifies discrepancies |
| **Pros/Cons** | 10 | Lists advantages and disadvantages of each method |
| **Recommendation** | 10 | Justified recommendation for which to use when |
| **Code Quality** | 10 | (Standard rubric, 10pt version) |
| **Visualization** | 5 | Comparison plots (e.g., scatter plot of counts) |

---

#### Assignment 3: QC Report

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **QC Metrics Calculation** | 15 | Calculates all essential metrics (genes, UMIs, mt%, etc.) |
| **Threshold Selection** | 20 | Clearly defines thresholds; provides rationale; uses appropriate method (MAD/fixed/quantile) |
| **Before/After Analysis** | 15 | Shows statistics and plots before and after filtering |
| **Doublet Detection** | 15 | Runs doublet detection; interprets scores; removes doublets |
| **Interpretation** | 15 | Assesses data quality; validates known markers; checks for over-filtering |
| **Code Quality** | 10 | (Standard rubric, 10pt version) |
| **Visualization** | 10 | QC violin plots, scatter plots, before/after comparisons |

---

### Course 2: Statistics & Differential Expression

#### Assignment 1: Design Matrix Construction

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Scenario 1** | 20 | Correct design matrix; contrast definition; explanation |
| **Scenario 2** | 20 | Correct design matrix; contrast definition; explanation |
| **Scenario 3** | 20 | Correct design matrix with interaction; contrast definition; explanation |
| **Confounding Discussion** | 15 | Identifies potential confounders; proposes solutions or caveats |
| **Code** | 15 | Shows how to construct matrices in R/Python |
| **Interpretation** | 10 | Explains what each coefficient represents |

---

#### Assignment 2: Multiple Testing Simulation

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Simulation Design** | 15 | Correctly simulates null data; appropriate parameters |
| **Raw p-value Analysis** | 15 | Counts false positives without correction; shows inflation |
| **FDR Correction** | 20 | Applies Benjamini-Hochberg; shows reduced false positives |
| **Threshold Comparison** | 15 | Compares p<0.05 vs FDR<0.05 vs FDR<0.1 |
| **Visualization** | 15 | Histogram of p-values; plot of # discoveries vs FDR threshold |
| **Code Quality** | 10 | (Standard rubric, 10pt version) |
| **Interpretation** | 10 | Explains multiple testing problem; recommends FDR threshold |

---

#### Assignment 3: Complete DESeq2 Analysis

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Data Loading** | 5 | Correctly loads count matrix and metadata |
| **DESeq2 Setup** | 10 | Proper design formula; creates DESeqDataSet |
| **DE Analysis** | 15 | Runs DESeq2 correctly; extracts results |
| **Results Filtering** | 10 | Applies FDR and log2FC thresholds; justifies choices |
| **Visualization** | 20 | MA plot, volcano plot, heatmap of top DEGs (all well-labeled) |
| **Gene Annotation** | 10 | Adds gene symbols/names to results |
| **Interpretation** | 15 | Biological interpretation of top DEGs; connects to known biology |
| **Code Quality** | 10 | (Standard rubric, 10pt version) |
| **Reproducibility** | 5 | Complete workflow; session info provided |

---

#### Assignment 4: Pseudobulk vs Cell-Level Comparison

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Pseudobulk Aggregation** | 15 | Correctly aggregates cells by sample; creates count matrix |
| **Pseudobulk DE** | 15 | Runs DESeq2 on pseudobulk data; proper design |
| **Cell-Level DE** | 15 | Runs Wilcoxon or equivalent on single-cell data |
| **Comparison** | 20 | Compares results; calculates overlap; shows discrepancies |
| **P-value Analysis** | 10 | Compares p-value distributions between methods |
| **Discussion** | 15 | Explains why methods differ; recommends appropriate approach |
| **Visualization** | 10 | Venn diagram, scatter plot of log2FC, p-value histograms |

---

### Course 3: Clustering & Annotation

#### Assignment 1: Dimensionality Reduction Comparison

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **PCA** | 15 | Runs PCA; selects appropriate # of PCs; visualizes variance explained |
| **UMAP** | 15 | Generates UMAP with ≥2 parameter sets; compares results |
| **t-SNE** | 15 | Generates t-SNE with ≥2 perplexity values; compares results |
| **Comparison** | 20 | Systematically compares methods; discusses trade-offs |
| **Visualization** | 20 | Side-by-side embeddings; elbow plot; parameter comparison |
| **Code Quality** | 10 | (Standard rubric, 10pt version) |
| **Interpretation** | 5 | Recommends method and parameters; justifies choice |

---

#### Assignment 2: Clustering at Multiple Resolutions

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Graph Construction** | 10 | Builds nearest neighbor graph correctly |
| **Multi-Resolution Clustering** | 20 | Clusters at 3-5 different resolutions; systematic approach |
| **Cluster Evaluation** | 20 | Calculates silhouette scores or other metrics at each resolution |
| **Marker Analysis** | 20 | Identifies markers for clusters at chosen resolution |
| **Resolution Selection** | 15 | Justifies chosen resolution based on markers and metrics |
| **Visualization** | 10 | UMAPs at different resolutions; marker heatmap/dotplot |
| **Interpretation** | 5 | Discusses trade-offs; identifies over/under-clustering |

---

#### Assignment 3: Complete Cell Type Annotation

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Marker Identification** | 15 | Finds markers for all clusters; appropriate statistical test |
| **Manual Annotation** | 20 | Assigns labels based on markers; uses canonical markers |
| **Automated Annotation** | 15 | Runs ≥1 automated tool (SingleR, CellTypist, etc.) |
| **Comparison** | 10 | Compares manual vs automated; resolves disagreements |
| **Confidence Assessment** | 15 | Assigns confidence levels; provides evidence table |
| **Visualization** | 15 | Annotated UMAP; marker dotplot; confidence indicators |
| **Documentation** | 10 | Evidence for each annotation; alternative interpretations noted |

---

### Course 4: Data Integration

#### Assignment 1: Batch Effect Diagnosis

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Visual Diagnosis** | 20 | UMAPs showing batch separation; PCA plots colored by batch |
| **Quantitative Metrics** | 20 | Calculates LISI, kBET, or similar metrics before integration |
| **Per-Cell-Type Analysis** | 15 | Checks batch effects within each cell type separately |
| **Severity Assessment** | 15 | Classifies batch effect as mild/moderate/severe; justifies |
| **Recommendation** | 15 | Recommends whether integration is needed; which method |
| **Visualization** | 10 | Clear plots showing batch effects |
| **Interpretation** | 5 | Discusses biological vs technical variation |

---

#### Assignment 2: Method Comparison

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Method 1 (Harmony)** | 15 | Runs Harmony correctly; appropriate parameters |
| **Method 2 (Seurat)** | 15 | Runs Seurat integration correctly |
| **Mixing Metrics** | 20 | Calculates batch mixing metrics for both methods |
| **Conservation Metrics** | 20 | Calculates biological conservation metrics for both |
| **Comparison Table** | 10 | Side-by-side comparison of methods and metrics |
| **Visualization** | 10 | UMAPs for both methods; metrics comparison plots |
| **Recommendation** | 10 | Recommends method based on quantitative evaluation |

---

#### Assignment 3: Full Integration Evaluation Report

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Before Integration** | 10 | Documents batch effects pre-integration |
| **Integration Execution** | 15 | Runs chosen method with justified parameters |
| **After Integration** | 15 | Shows improved batch mixing post-integration |
| **Quantitative Evaluation** | 25 | Comprehensive metrics (mixing + conservation) |
| **Marker Preservation** | 15 | Validates that known markers still identify cell types |
| **Visualization** | 10 | Before/after UMAPs; metrics plots; marker plots |
| **Discussion** | 10 | Assesses integration quality; identifies any over-correction |

---

### Course 5: Trajectory Analysis

#### Assignment 1: Pseudotime Inference Comparison

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Method 1 (DPT/PAGA)** | 15 | Runs diffusion pseudotime correctly |
| **Method 2 (Monocle/Slingshot)** | 15 | Runs alternative method correctly |
| **Root Selection** | 15 | Justifies root cell selection; tests alternatives if ambiguous |
| **Comparison** | 20 | Compares pseudotime orderings; calculates correlation |
| **Marker Validation** | 15 | Checks that known temporal markers change as expected |
| **Visualization** | 15 | Pseudotime on UMAP; marker trends along pseudotime |
| **Interpretation** | 5 | Discusses agreement/disagreement; recommends approach |

---

#### Assignment 2: RNA Velocity Analysis

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Velocity Calculation** | 20 | Runs scVelo correctly; appropriate model (stochastic/dynamical) |
| **Quality Assessment** | 15 | Checks velocity confidence; gene-level diagnostics |
| **Visualization** | 20 | Velocity stream plot; velocity embedding; gene plots |
| **Direction Validation** | 15 | Validates direction with known stage markers |
| **Interpretation** | 15 | Interprets velocity patterns; notes limitations |
| **Code Quality** | 10 | (Standard rubric, 10pt version) |
| **Discussion** | 5 | Discusses when velocity is reliable vs uncertain |

---

#### Assignment 3: Trajectory DE and Visualization

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Trajectory DE** | 25 | Identifies genes changing along pseudotime; appropriate method |
| **Gene Clustering** | 15 | Groups genes by expression pattern (e.g., early/late) |
| **Visualization** | 25 | Expression heatmap along pseudotime; individual gene trends; module plots |
| **Biological Interpretation** | 20 | Relates gene patterns to known biology; pathway enrichment |
| **Code Quality** | 10 | (Standard rubric, 10pt version) |
| **Documentation** | 5 | Clear workflow; parameter justification |

---

### Course 6: Cell-Cell Communication

#### Assignment 1: Data Preparation for CCC

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Annotation Quality Check** | 20 | Validates cell type annotations; sufficient confidence |
| **Expression Thresholds** | 15 | Sets appropriate thresholds for L-R detection |
| **Sample Structure** | 15 | Documents sample/batch structure; notes confounders |
| **Data Formatting** | 20 | Prepares data in correct format for CCC tools |
| **Positive Controls** | 15 | Identifies known interactions that should be detected |
| **Code Quality** | 10 | (Standard rubric, 10pt version) |
| **Documentation** | 5 | Clear rationale for all decisions |

---

#### Assignment 2: Tool Comparison

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **CellPhoneDB** | 15 | Runs CellPhoneDB correctly; interprets output |
| **CellChat** | 15 | Runs CellChat correctly; extracts key pathways |
| **Comparison Analysis** | 25 | Compares top interactions from both tools; calculates overlap |
| **Agreement Analysis** | 15 | Identifies consistently detected vs tool-specific interactions |
| **Visualization** | 15 | Dotplots, network plots, Venn diagrams |
| **Code Quality** | 10 | (Standard rubric, 10pt version) |
| **Discussion** | 5 | Explains why tools agree/disagree; recommends approach |

---

#### Assignment 3: NicheNet + LIANA Consensus

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **NicheNet Analysis** | 20 | Runs NicheNet; identifies ligand activities |
| **LIANA Consensus** | 20 | Runs LIANA with multiple methods; ranks interactions |
| **Integration** | 15 | Integrates NicheNet ligand-target with LIANA interactions |
| **Prioritization** | 20 | Creates prioritized list of interactions for validation |
| **Visualization** | 15 | Ligand activity plot; consensus ranking; target gene networks |
| **Code Quality** | 5 | (Standard rubric, 5pt version) |
| **Interpretation** | 5 | Biological interpretation; validation suggestions |

---

#### Assignment 4: CCC Network Analysis

**Total: 100 points**

| Component | Points | Criteria |
|-----------|--------|----------|
| **Network Construction** | 20 | Builds CCC network (nodes=cell types, edges=interactions) |
| **Centrality Analysis** | 20 | Calculates degree, betweenness, etc.; identifies hubs |
| **Community Detection** | 15 | Identifies modules/communities in CCC network |
| **Visualization** | 20 | Network diagram; centrality plots; community plots |
| **Biological Interpretation** | 15 | Interprets network structure in biological context |
| **Code Quality** | 5 | (Standard rubric, 5pt version) |
| **Discussion** | 5 | Discusses limitations; proposes validation experiments |

---

## Final Project Rubrics

### Final Project - General Template (100 points)

| Component | Points | Criteria |
|-----------|--------|----------|
| **Data Processing** | 15 | Complete workflow from raw to processed data |
| **Analysis Execution** | 25 | All required analyses performed correctly |
| **Integration** | 15 | Integrates concepts from multiple course modules |
| **Interpretation** | 20 | Comprehensive biological interpretation |
| **Visualization** | 10 | Publication-quality figures throughout |
| **Code Quality** | 5 | Well-organized, documented code |
| **Reproducibility** | 5 | Complete repository with README; can be re-run |
| **Written Report** | 5 | Clear, concise 1-2 page biological interpretation |

---

## Grading Guidelines

### Late Submissions
- **1-24 hours late:** -10%
- **24-48 hours late:** -20%
- **48-72 hours late:** -30%
- **>72 hours late:** Case-by-case (contact instructor)

### Collaboration Policy
- **Encouraged:** Discussing concepts, helping debug
- **Not Allowed:** Sharing code directly, copying analyses
- **Required:** Cite any code adapted from external sources
- **Suspected plagiarism:** 0 points + academic integrity review

### Resubmission Policy
- **First submission:** Full rubric applies
- **Resubmission (if allowed):** Maximum 85% of points
- **Must clearly mark changes** in resubmission

---

## Example Grading Comments

### Excellent Work (90-100)
> "Excellent work! Your code is well-organized and clearly documented. The interpretation goes beyond the immediate results to connect with broader biological context. The visualizations are publication-quality. Minor suggestion: consider adding a sensitivity analysis for the threshold choice."

### Good Work (80-89)
> "Good job overall. The analysis is correct and well-executed. Documentation could be improved - add more comments explaining why you chose specific parameters. The interpretation is solid but could benefit from discussing limitations. Figures are clear but labels could be more descriptive."

### Satisfactory Work (70-79)
> "Satisfactory work. The basic analyses are correct, but there are several areas for improvement: 1) Code organization could be cleaner with better function usage, 2) More documentation needed - why these thresholds?, 3) Interpretation is present but lacks depth, 4) Some figures are missing axis labels. Please review the rubric for specific improvement areas."

### Needs Improvement (<70)
> "This submission needs significant revision. Major issues: 1) Code has errors that prevent full execution, 2) Key analyses are missing or incorrect, 3) Little to no interpretation provided, 4) Figures are unclear or mislabeled. Please review course materials and office hours before resubmission."

---

## Self-Assessment Checklist

Before submitting, ask yourself:

### Code Quality
- [ ] Can someone else understand my code?
- [ ] Are variables named descriptively?
- [ ] Is the code organized logically?
- [ ] Does it run without errors?

### Documentation
- [ ] Have I explained my reasoning?
- [ ] Are parameters justified?
- [ ] Can someone reproduce this?
- [ ] Are versions recorded?

### Correctness
- [ ] Did I use appropriate methods?
- [ ] Did I validate my results?
- [ ] Are there any obvious errors?
- [ ] Do the results make sense?

### Interpretation
- [ ] Have I explained results in biological terms?
- [ ] Have I acknowledged limitations?
- [ ] Are conclusions supported by results?
- [ ] Have I considered alternatives?

### Visualization
- [ ] Are all axes labeled?
- [ ] Are legends present and clear?
- [ ] Are colors appropriate?
- [ ] Do figures support my conclusions?

---

## Additional Resources

### Code Style Guides
- **Python:** [PEP 8](https://pep8.org/)
- **R:** [Tidyverse Style Guide](https://style.tidyverse.org/)

### Example Assignments
- Each course folder contains `/examples/` with annotated example submissions

### Getting Help
- Office hours: [Schedule]
- Discussion forum: [Link]
- Peer review groups: [Sign up]

---

**Last Updated:** February 2, 2026  
**Questions?** Contact course instructors or post in discussion forum
