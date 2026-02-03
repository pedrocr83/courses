# Self-Assessment Checklist: Course 1 - Single-Cell Processing

**Track your mastery of key concepts and skills**

---

## How to Use This Checklist

- Review after each module
- Check off items you can confidently do/explain
- Revisit unchecked items before moving forward
- Use as study guide for quizzes and assignments

**Mastery Levels:**
- ✅ = Can do confidently
- 🔄 = Need more practice
- ❌ = Don't understand yet

---

## Module 1: Introduction to scRNA-seq Processing

### Concepts
- [ ] I can explain why processing quality determines analysis quality
- [ ] I can name the main stages of an scRNA-seq pipeline
- [ ] I can describe differences between 10x, Smart-seq2, and Drop-seq
- [ ] I understand what a count matrix is and its structure

### Skills
- [ ] I can identify the components of a published scRNA-seq dataset
- [ ] I can navigate a typical scRNA-seq project directory structure

### Red Flags I Can Spot
- [ ] Missing metadata
- [ ] Unorganized file structure
- [ ] Unknown platform/chemistry

---

## Module 2: Experimental Design & Metadata

### Concepts
- [ ] I can distinguish biological from technical replicates
- [ ] I understand how batch effects arise
- [ ] I know what metadata to track and why
- [ ] I can identify confounded experimental designs

### Skills
- [ ] I can create a proper metadata template for an experiment
- [ ] I can identify critical missing metadata
- [ ] I can plan a batch-balanced experimental design

### Red Flags I Can Spot
- [ ] Batch confounded with condition
- [ ] No replicates
- [ ] Missing sample information

---

## Module 3: FASTQ Files & Read Structure

### Concepts
- [ ] I understand FASTQ file format
- [ ] I can identify cell barcode, UMI, and cDNA in reads
- [ ] I know how quality scores work
- [ ] I understand protocol-specific read configurations

### Skills
- [ ] I can inspect FASTQ files from command line
- [ ] I can identify 10x v2 vs v3 chemistry from read structure
- [ ] I can check FASTQ quality
- [ ] I can estimate expected cell count from FASTQ

### Red Flags I Can Spot
- [ ] Very low quality scores
- [ ] Unexpected read length
- [ ] Barcode structure doesn't match protocol

---

## Module 4: Reference Preparation

### Concepts
- [ ] I understand what a reference genome is
- [ ] I know what a GTF annotation file contains
- [ ] I understand version matching importance
- [ ] I know when custom references are needed

### Skills
- [ ] I can download reference files from Ensembl/GENCODE
- [ ] I can build a Cell Ranger reference
- [ ] I can verify reference integrity
- [ ] I can document reference versions

### Red Flags I Can Spot
- [ ] Wrong species reference
- [ ] Mismatched FASTA and GTF versions
- [ ] Missing genes in annotation

---

## Module 5: Alignment & Quantification

### Concepts
- [ ] I can explain full alignment vs pseudo-alignment
- [ ] I understand the Cell Ranger workflow
- [ ] I know alternatives to Cell Ranger (STARsolo, kb)
- [ ] I understand mapping rate interpretation

### Skills
- [ ] I can run Cell Ranger count
- [ ] I can run kallisto|bustools
- [ ] I can interpret pipeline metrics
- [ ] I can compare pipeline outputs

### Red Flags I Can Spot
- [ ] Very low mapping rate (<50%)
- [ ] High multi-mapping percentage
- [ ] Unexpectedly low gene detection

---

## Module 6: Cell Calling & UMI Deduplication

### Concepts
- [ ] I understand the cell calling problem
- [ ] I can explain knee plot interpretation
- [ ] I understand EmptyDrops method
- [ ] I know how UMI collapsing works

### Skills
- [ ] I can generate and interpret a knee plot
- [ ] I can apply different cell calling methods
- [ ] I can identify good vs poor cell calling
- [ ] I can troubleshoot cell calling failures

### Red Flags I Can Spot
- [ ] No clear inflection point
- [ ] Cell count wildly different from expectation
- [ ] Many "cells" with very low UMI counts

---

## Module 7: Count Matrix

### Concepts
- [ ] I understand count matrix dimensions (genes × cells)
- [ ] I know what sparse matrices are and why they're used
- [ ] I can distinguish raw vs filtered matrices
- [ ] I understand matrix file formats (MTX, H5, H5AD)

### Skills
- [ ] I can load 10x MTX format in Python/R
- [ ] I can load H5AD files
- [ ] I can check matrix sparsity
- [ ] I can extract subsets of data

### Red Flags I Can Spot
- [ ] Unexpected matrix dimensions
- [ ] Missing barcodes or features files
- [ ] Corrupted H5 files

---

## Module 8: Quality Control Metrics

### Concepts
- [ ] I understand key QC metrics (n_genes, total_counts, pct_mt)
- [ ] I know what high mt% indicates
- [ ] I can interpret QC distributions
- [ ] I understand tissue-specific QC considerations

### Skills
- [ ] I can calculate all essential QC metrics
- [ ] I can create violin plots and scatter plots
- [ ] I can identify QC outliers
- [ ] I can apply the QC Decision Framework

### Red Flags I Can Spot
- [ ] Bimodal distributions without biological reason
- [ ] Very high mt% (>25% in non-neurons)
- [ ] Very low gene detection across samples

---

## Module 9: Filtering Strategies

### Concepts
- [ ] I understand MAD-based thresholding
- [ ] I know the dangers of over-filtering
- [ ] I understand doublet detection principles
- [ ] I can use Scrublet or DoubletFinder

### Skills
- [ ] I can set QC thresholds using MAD method
- [ ] I can apply tissue-specific thresholds
- [ ] I can run doublet detection
- [ ] I can validate filtering decisions

### Red Flags I Can Spot
- [ ] Filtering removes >50% of cells
- [ ] Expected cell types disappear
- [ ] QC metrics still drive clustering after filtering

---

## Module 10: Data Export & Reproducibility

### Concepts
- [ ] I understand different output formats (MTX, H5, H5AD, Loom)
- [ ] I know what to save and never overwrite
- [ ] I understand reproducibility requirements
- [ ] I can document processing pipelines

### Skills
- [ ] I can export data in multiple formats
- [ ] I can create a processing log
- [ ] I can record all software versions
- [ ] I can make my analysis reproducible

### Red Flags I Can Spot
- [ ] Overwriting raw data
- [ ] Missing version information
- [ ] Undocumented parameter choices

---

## Overall Competency Assessment

### Basic Level (Can proceed to Course 2)
- [ ] I can run a complete processing pipeline
- [ ] I can perform basic QC
- [ ] I understand major concepts
- [ ] I can troubleshoot common errors
- [ ] **Score:** ≥70% of items checked

### Proficient Level (Ready for independent work)
- [ ] I can choose appropriate methods for different scenarios
- [ ] I can troubleshoot complex problems
- [ ] I can justify all processing decisions
- [ ] I can document reproducibly
- [ ] **Score:** ≥85% of items checked

### Expert Level (Can teach others)
- [ ] I understand statistical foundations
- [ ] I can optimize workflows
- [ ] I can handle edge cases
- [ ] I can critically evaluate methods
- [ ] **Score:** ≥95% of items checked

---

## Action Items for Unchecked Items

### If you have many unchecked items:
1. Review relevant module in `course_single_cell_processing.md`
2. Redo the corresponding lab
3. Watch supplementary videos (if available)
4. Practice with different datasets
5. Ask for help in discussion forum

### Priority order for re-study:
1. **Critical:** Modules 6, 8, 9 (cell calling, QC, filtering)
2. **Important:** Modules 5, 7, 10 (alignment, matrices, export)
3. **Foundational:** Modules 1-4 (background concepts)

---

## Before Each Assignment

### Assignment 1 Prerequisites
- [ ] Modules 1-3 mostly checked
- [ ] Can inspect FASTQ files
- [ ] Can identify read structure

### Assignment 2 Prerequisites
- [ ] Module 5 mostly checked
- [ ] Can run both Cell Ranger and kb
- [ ] Can compare pipeline outputs

### Assignment 3 Prerequisites
- [ ] Modules 8-9 fully checked
- [ ] Can apply QC Decision Framework
- [ ] Can run and interpret doublet detection

### Final Project Prerequisites
- [ ] ALL modules ≥80% checked
- [ ] Can execute end-to-end pipeline
- [ ] Can document reproducibly

---

## Study Tips

### For Concept Understanding:
- Draw flowcharts of pipelines
- Explain concepts to someone else
- Create comparison tables
- Read original method papers

### For Skill Development:
- Practice on multiple datasets
- Try breaking things and fixing them
- Vary parameters and observe effects
- Create your own examples

### For Troubleshooting:
- Keep a log of errors encountered
- Document solutions that worked
- Build your own troubleshooting guide
- Help others (teaching reinforces learning)

---

## Progress Tracking

**Date started:** _______________  
**Target completion:** _______________

**Weekly check-ins:**
- Week 1: ____% complete
- Week 2: ____% complete
- Week 3: ____% complete

**Areas needing extra focus:**
1. _______________
2. _______________
3. _______________

**Questions for instructor:**
1. _______________
2. _______________
3. _______________

---

## Course Completion Checklist

Before claiming course completion:
- [ ] All self-assessment items ≥80% checked
- [ ] All labs completed
- [ ] All assignments submitted
- [ ] All quizzes passed (≥70%)
- [ ] Final project completed
- [ ] Can confidently process a new dataset independently

**Ready for Course 2?** ✅

---

**Remember:** This is for YOUR learning. Be honest with yourself—unchecked items are opportunities for growth, not failures!
