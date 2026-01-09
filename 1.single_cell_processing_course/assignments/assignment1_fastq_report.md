# Assignment 1: FASTQ Inspection Report

**Module:** 3 - FASTQ Files & Read Structure  
**Due:** End of Week 1  
**Points:** 100

---

## Overview

In this assignment, you will inspect FASTQ files from a 10x Genomics scRNA-seq experiment and write a short report documenting your findings.

---

## Tasks

### Part A: FASTQ Structure (40 points)

1. Download the 1k PBMC demo FASTQ files (see setup.md)
2. Inspect the first 10 reads from both R1 and R2 files
3. Answer the following:
   - What is the length of R1 reads?
   - What is the length of R2 reads?
   - Identify the cell barcode and UMI positions in R1

### Part B: Quality Scores (30 points)

1. Extract quality scores from 100 reads
2. Calculate the average quality score per position
3. Create a simple plot showing quality across read positions
4. Identify if there's a quality drop-off pattern

### Part C: Read Count Summary (30 points)

1. Count total reads in one FASTQ file (use `zcat file.fastq.gz | wc -l` divided by 4)
2. Calculate the approximate sequencing depth per cell (assuming 1,000 cells)
3. Discuss: Is this sufficient depth for typical scRNA-seq analysis?

---

## Deliverables

1. **Report (PDF or Markdown):**
   - 1-2 pages
   - Include code snippets used
   - Include quality plot

2. **Code:**
   - Python script or Jupyter notebook with your analysis

---

## Rubric

| Criterion | Points |
|-----------|--------|
| Correct FASTQ structure identification | 20 |
| Accurate barcode/UMI positions | 20 |
| Quality score analysis | 20 |
| Quality plot | 10 |
| Read count and depth calculation | 15 |
| Discussion of depth adequacy | 10 |
| Code quality and documentation | 5 |

---

## Submission

Submit your report and code to the course assignments folder.

