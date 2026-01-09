# Quiz: Module 1 - Introduction to scRNA-seq Processing

**10 Questions | Passing: 70%**

---

## Question 1
What is the final output of the scRNA-seq processing pipeline?

- A) FASTQ files
- B) BAM alignment files
- C) Cell × Gene count matrix
- D) Cluster assignments

---

## Question 2
Which statement best describes why processing quality matters?

- A) It only affects storage requirements
- B) It sets the ceiling for all downstream analysis quality
- C) It determines the number of cells captured
- D) It only matters for publication

---

## Question 3
What does "UMI" stand for?

- A) Universal Mapping Index
- B) Unique Molecular Identifier
- C) Unified Matrix Input
- D) Unit of Measurement Index

---

## Question 4
Which platform is most commonly used for droplet-based scRNA-seq?

- A) Illumina HiSeq
- B) Smart-seq2
- C) 10x Genomics
- D) PacBio

---

## Question 5
What is the primary purpose of a cell barcode?

- A) Identify which gene a read came from
- B) Identify which cell a read came from
- C) Measure sequencing quality
- D) Remove PCR duplicates

---

## Question 6
In a typical scRNA-seq experiment, which is larger?

- A) Number of genes
- B) Number of cells
- C) They are always equal
- D) Depends entirely on the tissue

---

## Question 7
What type of noise does scRNA-seq processing aim to minimize?

- A) Only biological variation
- B) Only technical artifacts
- C) Both biological variation and technical artifacts
- D) Neither - noise cannot be minimized

---

## Question 8
What happens to data quality if empty droplets are incorrectly called as cells?

- A) No effect
- B) Artificially inflates cell count with low-quality "cells"
- C) Reduces the number of genes detected
- D) Improves clustering results

---

## Question 9
Which of these is NOT a typical step in scRNA-seq processing?

- A) Alignment
- B) Cell clustering
- C) UMI deduplication
- D) Cell barcode detection

---

## Question 10
Why is metadata tracking critical in scRNA-seq experiments?

- A) It's only needed for publication
- B) Without it, data cannot be properly analyzed or reproduced
- C) It only matters for multi-sample experiments
- D) Metadata is optional

---

## Answer Key

1. C
2. B
3. B
4. C
5. B
6. A (typically ~20,000-30,000 genes vs hundreds to thousands of cells)
7. B (processing focuses on technical artifacts; biological variation is preserved)
8. B
9. B (clustering is downstream analysis, not processing)
10. B

