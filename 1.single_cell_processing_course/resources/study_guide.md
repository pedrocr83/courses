# Single-Cell RNA-seq Sample Processing

## From Raw Sequencing Reads to a High-Quality Expression Matrix

This document is a **practical, end-to-end study guide** for processing single-cell RNA-seq (scRNA-seq) data. It focuses on *what happens before analysis*: alignment, quantification, QC, filtering, and best practices. It is written for beginners but follows industry-grade standards.

---

## 1. What Is Single-Cell RNA-seq Processing?

Single-cell processing converts **raw sequencing reads (FASTQ)** into a **cell × gene count matrix** that accurately represents biological signal while minimizing technical noise.

The quality of everything downstream (DE, clustering, CCC, ML) depends on this step.

---

## 2. Overview of the scRNA-seq Pipeline

### High-Level Steps
1. Experimental design & metadata
2. Raw data (FASTQ)
3. Read alignment or pseudo-alignment
4. UMI collapsing & quantification
5. Cell barcode detection
6. Quality control (QC)
7. Filtering & cleanup
8. Final count matrix

---

## 3. Experimental Design & Metadata (Before Sequencing)

### Critical Decisions
- Platform (10x Genomics, Smart-seq2, Drop-seq, etc.)
- Number of cells vs sequencing depth
- Biological replicates
- Batch structure

### Metadata to Track
- Sample ID
- Donor
- Condition
- Batch / lane
- Chemistry version

Bad metadata = unusable data.

---

## 4. Raw Data: FASTQ Files

### What FASTQ Contains
- Read sequence
- Base quality scores

### scRNA-seq Specific Reads
- **Cell barcode** (which cell)
- **UMI** (which molecule)
- **cDNA read** (which gene)

Read structure depends on protocol.

---

## 5. Alignment vs Pseudo-Alignment

### Full Alignment
Maps reads to a reference genome.
- Accurate
- Slower

### Pseudo-Alignment
Maps reads to transcripts without base-level alignment.
- Faster
- Standard for scRNA-seq

Most modern pipelines use pseudo-alignment.

---

## 6. Reference Preparation

### Reference Components
- Genome (FASTA)
- Annotation (GTF)
- Transcriptome

### Best Practices
- Use matching genome + annotation versions
- Include spike-ins if applicable
- Track reference version in metadata

---

## 7. Popular Processing Pipelines

### Cell Ranger (10x Genomics)
- Industry standard
- End-to-end solution
- Highly reproducible

Docs: https://www.10xgenomics.com/support/software/cell-ranger

### STARsolo
- Based on STAR aligner
- Flexible, open-source

Docs: https://github.com/alexdobin/STAR

### kallisto | bustools
- Fast pseudo-alignment
- Modular design

Docs: https://www.kallistobus.tools/

---

## 8. Cell Barcode Detection

### The Problem
Which barcodes correspond to real cells vs empty droplets?

### Common Approaches
- Knee / inflection point
- EmptyDrops (statistical test)

Correct cell calling is critical.

---

## 9. UMI Deduplication

### Why UMIs Exist
- PCR duplicates inflate counts

### Deduplication
- Collapse reads with same (cell, gene, UMI)
- Correct sequencing errors in UMIs

This step defines molecule counts.

---

## 10. The Count Matrix

### Structure
- Rows: genes
- Columns: cells
- Values: UMI counts

This matrix is the *input* for all downstream analyses.

---

## 11. Quality Control (QC): Core Metrics

### Per-Cell Metrics
- Number of detected genes
- Total UMI counts
- % mitochondrial reads

### Per-Gene Metrics
- Number of cells detected

QC metrics detect dead cells, debris, and technical artifacts.

---

## 12. Mitochondrial Content

### Why It Matters
- High mitochondrial RNA → dying or stressed cells

### Typical Thresholds
- 5–20% depending on tissue

Thresholds are dataset-dependent.

---

## 13. Filtering Cells

### Common Filters
- Min genes per cell
- Max genes per cell (doublets)
- Max mitochondrial percentage

Filtering removes non-biological signal.

---

## 14. Doublet Detection

### What Are Doublets?
Two cells captured as one.

### Tools
- Scrublet
- DoubletFinder

Doublets distort clustering and DE.

---

## 15. Batch Effects (Early Awareness)

Batch effects originate during:
- Library prep
- Sequencing
- Sample handling

Good processing minimizes—but does not eliminate—batch effects.

---

## 16. Data Export & Formats

### Common Formats
- Matrix Market (.mtx)
- HDF5 (.h5, .h5ad)
- Loom

### Best Practices
- Store raw counts
- Store filtered counts
- Never overwrite raw data

---

## 17. Reproducibility Best Practices

- Fix software versions
- Log parameters
- Keep raw FASTQs
- Track pipeline configs

Processing must be reproducible.

---

## 18. Common Mistakes

- Using wrong reference genome
- Over-filtering cells
- Ignoring doublets
- Losing metadata

Most analysis problems start here.

---

## 19. Learning Path

### Beginner
- Run Cell Ranger on demo data
- Inspect QC metrics visually

### Intermediate
- Compare pipelines
- Tune QC thresholds

### Advanced
- Custom references
- Multi-modal data
- Large-scale automation

---

## 20. How This Connects to Downstream Analysis

Good processing enables:
- Reliable clustering
- Accurate differential expression
- Meaningful cell–cell communication
- Trustworthy ML models

Processing quality sets the ceiling for insight.

---

## 21. Suggested Outcome

By mastering this stage, you will:
- Understand where noise comes from
- Trust your count matrix
- Know when data is unusable
- Build robust downstream analyses

---

*This document is intended as the practical backbone of any single-cell transcriptomics curriculum and pairs naturally with statistics, DE, and cell–cell communication modules.*