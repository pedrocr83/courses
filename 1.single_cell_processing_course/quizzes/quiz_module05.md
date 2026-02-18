# Quiz: Module 5 - Alignment & Quantification Pipelines

**15 Questions | Passing: 70% | Mix: MC + Short Answer**

---

## Part A: Multiple Choice (10 questions)

### Question 1
What is the main difference between alignment and pseudo-alignment?

- A) Alignment is always faster
- B) Pseudo-alignment doesn't produce base-level mapping positions
- C) Alignment requires less memory
- D) Pseudo-alignment only works with 10x data

**Your answer:** _____

---

### Question 2
Which pipeline is the industry standard for 10x Genomics data?

- A) STAR
- B) kallisto
- C) Cell Ranger
- D) HTSeq

**Your answer:** _____

---

### Question 3
What does STARsolo provide that basic STAR does not?

- A) Genome alignment
- B) scRNA-seq specific features (barcode/UMI handling)
- C) Faster processing
- D) Better quality scores

**Your answer:** _____

---

### Question 4
kallisto | bustools is known for:

- A) Slow but accurate alignment
- B) Fast pseudo-alignment
- C) Only working with human data
- D) Requiring a GPU

**Your answer:** _____

---

### Question 5
What is the primary input to an alignment pipeline?

- A) Count matrix
- B) FASTQ files and reference
- C) BAM files
- D) H5AD files

**Your answer:** _____

---

### Question 6
Why might you choose pseudo-alignment over traditional alignment?

- A) More accurate gene quantification
- B) Faster processing with similar accuracy for quantification
- C) Better splice detection
- D) Required for droplet-based data

**Your answer:** _____

---

### Question 7
What metric indicates how well reads mapped to the reference?

- A) UMI count
- B) Mapping rate / alignment rate
- C) Cell count
- D) Sparsity

**Your answer:** _____

---

### Question 8
Which reference files are needed for most scRNA-seq pipelines?

- A) Only the genome FASTA
- B) Only the GTF annotation
- C) Both genome FASTA and GTF annotation
- D) Neither - references are built-in

**Your answer:** _____

---

### Question 9
What happens if you use the wrong genome assembly version?

- A) No effect on results
- B) Significant data loss from failed mapping
- C) Only affects visualization
- D) Pipeline automatically corrects it

**Your answer:** _____

---

### Question 10
Which is TRUE about Cell Ranger?

- A) It's open source and free
- B) It only works on Linux
- C) It provides an end-to-end solution with web reports
- D) It requires a GPU

**Your answer:** _____

---

## Part B: Short Answer (5 questions)

### Question 11
List the main reference inputs needed to run Cell Ranger or kallisto|bustools (be specific about file types).

**Your answer:**

---

### Question 12
In one sentence, explain why pseudo-alignment can be faster than traditional alignment while still producing accurate counts.

**Your answer:**

---

### Question 13
What is the purpose of the transcriptome index in kallisto, and how does it differ from a genome index?

**Your answer:**

---

### Question 14
Name two key outputs (or metrics) you would check after running an alignment pipeline to validate success.

**Your answer:**

---

### Question 15
When might you choose kallisto|bustools over Cell Ranger for 10x data?

**Your answer:**

---

## Answer Key

### Part A
1. B  
2. C  
3. B  
4. B  
5. B  
6. B  
7. B  
8. C  
9. B  
10. C  

### Part B (sample answers)
11. Genome FASTA (.fa), GTF annotation (.gtf), and optionally a pre-built transcriptome or index. 10x requires chemistry-specific reference.  
12. Pseudo-alignment finds which transcripts a read is compatible with without computing exact base positions, reducing computation while preserving count accuracy.  
13. Transcriptome index is built from transcript sequences (cDNA); kallisto aligns to transcripts directly. Genome index aligns to the full genome; gene assignment comes later.  
14. Mapping rate, number of cells detected, median genes/cell, total UMI count, fraction of reads in cells.  
15. When you need open-source/free software, faster runtime, or custom workflows; when avoiding 10x licensing.
