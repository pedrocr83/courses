# Quiz: Module 5 - Alignment & Quantification Pipelines

**10 Questions | Passing: 70%**

---

## Question 1
What is the main difference between alignment and pseudo-alignment?

- A) Alignment is always faster
- B) Pseudo-alignment doesn't produce base-level mapping positions
- C) Alignment requires less memory
- D) Pseudo-alignment only works with 10x data

---

## Question 2
Which pipeline is the industry standard for 10x Genomics data?

- A) STAR
- B) kallisto
- C) Cell Ranger
- D) HTSeq

---

## Question 3
What does STARsolo provide that basic STAR does not?

- A) Genome alignment
- B) scRNA-seq specific features (barcode/UMI handling)
- C) Faster processing
- D) Better quality scores

---

## Question 4
kallisto | bustools is known for:

- A) Slow but accurate alignment
- B) Fast pseudo-alignment
- C) Only working with human data
- D) Requiring a GPU

---

## Question 5
What is the primary input to an alignment pipeline?

- A) Count matrix
- B) FASTQ files and reference
- C) BAM files
- D) H5AD files

---

## Question 6
Why might you choose pseudo-alignment over traditional alignment?

- A) More accurate gene quantification
- B) Faster processing with similar accuracy for quantification
- C) Better splice detection
- D) Required for droplet-based data

---

## Question 7
What metric indicates how well reads mapped to the reference?

- A) UMI count
- B) Mapping rate / alignment rate
- C) Cell count
- D) Sparsity

---

## Question 8
Which reference files are needed for most scRNA-seq pipelines?

- A) Only the genome FASTA
- B) Only the GTF annotation
- C) Both genome FASTA and GTF annotation
- D) Neither - references are built-in

---

## Question 9
What happens if you use the wrong genome assembly version?

- A) No effect on results
- B) Significant data loss from failed mapping
- C) Only affects visualization
- D) Pipeline automatically corrects it

---

## Question 10
Which is TRUE about Cell Ranger?

- A) It's open source and free
- B) It only works on Linux
- C) It provides an end-to-end solution with web reports
- D) It requires a GPU

---

## Answer Key

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

