# Quiz: Module 1 - Introduction to scRNA-seq Processing

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What is the final output of the scRNA-seq processing pipeline?

- A) FASTQ files
- B) BAM alignment files
- C) Cell × Gene count matrix
- D) Cluster assignments

**Your answer:** _____

---

### Question 2
Which statement best describes why processing quality matters?

- A) It only affects storage requirements
- B) It sets the ceiling for all downstream analysis quality
- C) It determines the number of cells captured
- D) It only matters for publication

**Your answer:** _____

---

### Question 3
What does "UMI" stand for and what problem does it solve?

- A) Universal Mapping Index – identifies genes
- B) Unique Molecular Identifier – corrects for PCR amplification bias
- C) Unified Matrix Input – standardizes file format
- D) Unit of Measurement Index – normalizes counts

**Your answer:** _____

---

### Question 4
Which platform is most commonly used for droplet-based scRNA-seq?

- A) Illumina HiSeq
- B) Smart-seq2
- C) 10x Genomics
- D) PacBio

**Your answer:** _____

---

### Question 5
In a typical 10x experiment, what is the order of: cell barcode, UMI, cDNA?

- A) UMI, barcode, cDNA
- B) Barcode, UMI, cDNA
- C) cDNA, barcode, UMI
- D) Barcode, cDNA, UMI

**Your answer:** _____

---

### Question 6
What happens to data quality if empty droplets are incorrectly called as cells?

- A) No effect
- B) Artificially inflates cell count with low-quality "cells"
- C) Reduces the number of genes detected
- D) Improves clustering results

**Your answer:** _____

---

### Question 7
Which of these is NOT a typical step in scRNA-seq processing?

- A) Alignment or pseudo-alignment
- B) Cell clustering
- C) UMI deduplication
- D) Cell barcode detection

**Your answer:** _____

---

### Question 8
Why is metadata tracking critical in scRNA-seq experiments?

- A) It's only needed for publication
- B) Without it, data cannot be properly analyzed or reproduced
- C) It only matters for multi-sample experiments
- D) Metadata is optional

**Your answer:** _____

---

### Question 9
What is the primary purpose of a cell barcode?

- A) Identify which gene a read came from
- B) Identify which cell a read came from
- C) Measure sequencing quality
- D) Remove PCR duplicates

**Your answer:** _____

---

### Question 10
In a typical scRNA-seq count matrix, which dimension is usually larger?

- A) Number of genes
- B) Number of cells
- C) They are always equal
- D) Depends entirely on the tissue

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
List the main stages of an scRNA-seq processing pipeline (in order), from raw reads to final matrix.

**Your answer:**

---

### Question 12
Explain in 1–2 sentences why processing quality "sets the ceiling" for downstream analysis.

**Your answer:**

---

### Question 13
Name two differences between 10x Genomics and Smart-seq2 platforms (e.g., throughput, read structure, cost).

**Your answer:**

---

### Question 14
What type of noise does scRNA-seq processing aim to minimize, and what type should it preserve?

**Your answer:**

---

### Question 15
Give one example of a "red flag" in experimental design or metadata that could make data unusable.

**Your answer:**

---

## Answer Key

### Part A
1. C  
2. B  
3. B  
4. C  
5. B  
6. B  
7. B  
8. B  
9. B  
10. A  

### Part B (sample answers)
11. Raw FASTQ → alignment/pseudo-alignment → UMI deduplication → cell barcode detection → QC → filtering → count matrix  
12. Garbage in, garbage out. If processing introduces artifacts or fails to remove technical noise, clustering, DE, and all downstream steps will be compromised.  
13. 10x: droplet-based, high throughput, 3' or 5' capture, shorter reads. Smart-seq2: plate-based, full-length, lower throughput, higher cost per cell.  
14. Minimize: technical noise (empty droplets, doublets, low-quality cells, batch effects). Preserve: biological variation (real cell-to-cell differences).  
15. Examples: batch confounded with condition; no replicates; missing sample/condition metadata; unknown chemistry version.
