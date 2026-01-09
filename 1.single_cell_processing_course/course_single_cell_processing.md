# Course: Single-Cell RNA-seq Sample Processing

## From Raw Sequencing Reads to a High-Quality Expression Matrix

---

## Course Overview

| | |
|---|---|
| **Duration** | 3 weeks (15-20 hours total) |
| **Format** | Self-paced with hands-on labs |
| **Level** | Beginner to Intermediate |
| **Prerequisites** | Basic command line, familiarity with genomics concepts |

## Start Here (Do This First)

- **Step-by-step checklist**: `START_HERE.md`
- **Glossary**: `glossary.md`
- **Environment setup**: `setup.md`

### Learning Objectives

By the end of this course, you will be able to:
1. Process raw scRNA-seq FASTQ files into a count matrix
2. Apply industry-standard QC metrics and filtering strategies
3. Detect and remove technical artifacts (doublets, empty droplets, dead cells)
4. Export data in standard formats for downstream analysis
5. Build reproducible processing pipelines

---

## Course Structure

### Week 1: Foundations & Raw Data

#### Module 1: Introduction to scRNA-seq Processing (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 1.1 | Why processing quality matters | Lecture |
| 1.2 | The scRNA-seq pipeline overview | Lecture |
| 1.3 | Platforms: 10x, Smart-seq2, Drop-seq | Reading |
| **Lab 1** | Explore a published scRNA-seq dataset structure | Hands-on |

**Key Concepts:** Pipeline stages, platform differences, quality ceiling principle

---

#### Module 2: Experimental Design & Metadata (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 2.1 | Critical design decisions | Lecture |
| 2.2 | Metadata: what to track and why | Lecture |
| 2.3 | Batch structure planning | Case study |
| **Lab 2** | Create a metadata template for a hypothetical experiment | Hands-on |

**Key Concepts:** Replicates, batch effects origin, metadata standards

---

#### Module 3: FASTQ Files & Read Structure (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 3.1 | Anatomy of a FASTQ file | Lecture |
| 3.2 | scRNA-seq read structure: barcode, UMI, cDNA | Lecture |
| 3.3 | Protocol-specific read configurations | Reading |
| **Lab 3** | Inspect FASTQ files, identify read components | Hands-on |

**Key Concepts:** Cell barcodes, UMIs, quality scores, read pairs

---

### Week 2: Alignment, Quantification & Cell Calling

#### Module 4: Reference Preparation (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 4.1 | Genome, annotation, and transcriptome | Lecture |
| 4.2 | Matching versions and sources | Lecture |
| 4.3 | Custom references and spike-ins | Reading |
| **Lab 4** | Download and prepare a reference for Cell Ranger | Hands-on |

**Key Concepts:** FASTA, GTF, version matching, reference indexing

---

#### Module 5: Alignment & Quantification Pipelines (3 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 5.1 | Full alignment vs pseudo-alignment | Lecture |
| 5.2 | Cell Ranger deep dive | Lecture + Demo |
| 5.3 | STARsolo and kallisto\|bustools | Lecture |
| 5.4 | Choosing the right pipeline | Discussion |
| **Lab 5A** | Run Cell Ranger on 10x demo data | Hands-on |
| **Lab 5B** | Run kallisto\|bustools on the same data | Hands-on |

**Key Concepts:** Mapping rate, speed vs accuracy, reproducibility

---

#### Module 6: Cell Barcode Detection & UMI Deduplication (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 6.1 | The cell calling problem | Lecture |
| 6.2 | Knee plots and EmptyDrops | Lecture |
| 6.3 | UMI collapsing and error correction | Lecture |
| **Lab 6** | Visualize barcode rank plot, compare calling methods | Hands-on |

**Key Concepts:** Empty droplets, ambient RNA, PCR duplicates, molecule counting

---

### Week 3: Quality Control, Filtering & Output

#### Module 7: The Count Matrix (1 hr)
| Lesson | Content | Format |
|--------|---------|--------|
| 7.1 | Matrix structure: genes × cells | Lecture |
| 7.2 | Sparse matrix formats | Lecture |
| **Lab 7** | Load and explore a count matrix in Python/R | Hands-on |

**Key Concepts:** Sparsity, matrix dimensions, raw vs filtered

---

#### Module 8: Quality Control Metrics (2.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 8.1 | Per-cell metrics: genes, UMIs, mito% | Lecture |
| 8.2 | Per-gene metrics | Lecture |
| 8.3 | Visualizing QC distributions | Lecture + Demo |
| 8.4 | Mitochondrial content interpretation | Case study |
| **Lab 8** | Generate QC plots, identify outliers | Hands-on |

**Key Concepts:** Violin plots, scatter plots, threshold selection

---

#### Module 9: Filtering Strategies (2 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 9.1 | Filtering cells: thresholds and rationale | Lecture |
| 9.2 | The dangers of over-filtering | Lecture |
| 9.3 | Doublet detection with Scrublet/DoubletFinder | Lecture |
| **Lab 9** | Apply filters, run doublet detection, compare before/after | Hands-on |

**Key Concepts:** MAD-based thresholds, doublet scores, filtering documentation

---

#### Module 10: Data Export & Reproducibility (1.5 hrs)
| Lesson | Content | Format |
|--------|---------|--------|
| 10.1 | Output formats: MTX, H5, H5AD, Loom | Lecture |
| 10.2 | What to save and never overwrite | Lecture |
| 10.3 | Reproducibility checklist | Lecture |
| **Lab 10** | Export filtered matrix, create processing log | Hands-on |

**Key Concepts:** Raw vs processed, version tracking, pipeline configs

---

## Assessments

### Quizzes (after each module)
- 5-10 questions covering key concepts
- Immediate feedback

### Practical Assignments

| Assignment | Module | Description |
|------------|--------|-------------|
| **A1** | 3 | FASTQ inspection report |
| **A2** | 5 | Pipeline comparison: Cell Ranger vs kallisto |
| **A3** | 8-9 | QC report with filtering decisions |

### Final Project
**Process a real scRNA-seq dataset end-to-end**
- Start from FASTQ
- Document all decisions
- Produce filtered count matrix
- Write 1-page QC summary

---

## Lab Environment Setup

### Required Software
```
# Core tools
- Cell Ranger 8.x
- STAR 2.7.x (for STARsolo)
- kallisto + bustools
- Python 3.10+ with scanpy
- R 4.x with Seurat (optional)

# Supporting tools
- samtools
- FastQC
```

### Demo Datasets
| Dataset | Source | Use |
|---------|--------|-----|
| 1k PBMCs | 10x Genomics | Labs 5-10 |
| 5k PBMCs | 10x Genomics | Final project option |

---

## Resources

### Primary References
- 10x Genomics Cell Ranger documentation
- Scanpy tutorials (Theis Lab)
- Orchestrating Single-Cell Analysis with Bioconductor

### Tools Documentation
| Tool | URL |
|------|-----|
| Cell Ranger | https://www.10xgenomics.com/support/software/cell-ranger |
| STARsolo | https://github.com/alexdobin/STAR |
| kallisto\|bustools | https://www.kallistobus.tools/ |
| Scrublet | https://github.com/swolock/scrublet |
| DoubletFinder | https://github.com/chris-mcginnis-ucsf/DoubletFinder |

---

## Common Mistakes to Avoid

| Mistake | Consequence | Prevention |
|---------|-------------|------------|
| Wrong reference genome | Massive data loss | Verify species + version |
| Over-filtering cells | Lose rare populations | Use adaptive thresholds |
| Ignoring doublets | Distorted clusters | Always run detection |
| Poor metadata | Unusable experiment | Template from day 1 |
| Overwriting raw data | Irreversible loss | Separate raw/processed dirs |

---

## What Bad Processing / Bad QC Looks Like (and How to Catch It)

### Symptoms (Red Flags)

- **FASTQ / read structure**: unexpected read lengths; barcode/UMI not where you think; barcode whitelist mismatch.
- **Alignment/quantification**: very low mapping rate; high multi-mapping; extreme rRNA/mt fraction; suspiciously low genes detected.
- **Cell calling**: no clear knee/inflection; “cells” dominated by low counts; cell count wildly different from expectations.
- **QC distributions**: bimodal `total_counts` without biological reason; very high `pct_counts_mt`; very low `n_genes_by_counts` for most cells.
- **Doublets**: clusters with mixed canonical markers (e.g., T + B markers), inflated UMIs/genes.

### Common Failure Modes → What It Looks Like → Fix

- **Wrong reference (species/version/annotation mismatch)**  
  - **Looks like**: low mapping, many unassigned reads, missing expected markers.  
  - **Fix**: rebuild reference with correct FASTA/GTF versions; record versions in `templates/processing_log_template.md`.

- **Ambient RNA / empty droplets called as cells**  
  - **Looks like**: lots of “cells” with very low UMIs/genes; markers appear weak everywhere.  
  - **Fix**: use knee plot + EmptyDrops-style logic; raise min counts; validate with marker enrichment.

- **Over-filtering**  
  - **Looks like**: rare populations disappear; strong shifts in cell-type proportions; QC thresholds too aggressive.  
  - **Fix**: choose thresholds per tissue; use MAD/outlier-based rules; justify thresholds in QC report.

- **Under-filtering**  
  - **Looks like**: QC plots show long low-quality tails; downstream clustering dominated by mt/ribo.  
  - **Fix**: remove dead/dying cells (high mito), remove very low gene cells, remove extreme high-UMI doublet candidates.

- **Barcode/chemistry mismatch (10x v2 vs v3 etc.)**  
  - **Looks like**: cell calling fails; UMI patterns odd; very low usable reads.  
  - **Fix**: confirm chemistry/read structure; ensure pipeline parameters match protocol.

### Minimum “Sanity Checks” Before Moving On

- **Counts**: cells and genes are plausible for the experiment.
- **Markers**: canonical markers exist in expected populations (even pre-clustering, spot check).
- **QC**: filtering decisions are documented and plots saved (before/after).
- **Reproducibility**: you can re-run from raw → matrix using recorded software versions and parameters.

---

## Learning Path Progression

```
BEGINNER (Weeks 1-3)
│
├── Run Cell Ranger on demo data
├── Inspect QC metrics visually
└── Apply standard filters
        │
        ▼
INTERMEDIATE (Post-course)
│
├── Compare multiple pipelines
├── Tune QC thresholds per tissue
├── Handle multi-sample experiments
└── Automate with Snakemake/Nextflow
        │
        ▼
ADVANCED (6+ months)
│
├── Build custom references
├── Process multi-modal data (CITE-seq, ATAC)
├── Large-scale automation (cloud)
└── Contribute to tool development
```

---

## Connection to Other Courses

| Prerequisite for | Why |
|------------------|-----|
| Differential Expression | Requires clean count matrix |
| Cell-Cell Communication | Cell identity depends on QC |
| Clustering & Annotation | Garbage in, garbage out |
| Machine Learning on scRNA | Feature quality matters |

---

## Suggested Schedule

| Day | Activity | Time |
|-----|----------|------|
| Mon/Wed/Fri | Lectures + Reading | 1 hr |
| Tue/Thu | Labs | 1.5 hrs |
| Weekend | Assignments | 2 hrs |

*Total: ~6 hours/week for 3 weeks*

---

## Certificate Criteria

To complete this course:
- [ ] Pass all module quizzes (≥70%)
- [ ] Submit all 3 practical assignments
- [ ] Complete final project with passing QC metrics
- [ ] Document reproducible pipeline

---

*This course prepares you for any downstream single-cell analysis and establishes the foundation for reliable, publishable results.*

