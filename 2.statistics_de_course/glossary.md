# Glossary: Statistics & Differential Expression

Use this while working through the Rmd labs and writing your reports.

---

## Core Statistical Terms

- **Count data**: Non-negative integers (reads/UMIs) measured per gene per sample/cell.
- **Poisson distribution**: Count model where mean = variance; often too simple for RNA-seq.
- **Overdispersion**: Variance > mean; common in RNA-seq counts.
- **Negative binomial (NB)**: Count model that handles overdispersion (core to DESeq2/edgeR).
- **GLM (Generalized Linear Model)**: Regression framework used for NB models in DESeq2/edgeR.
- **Dispersion**: Parameter controlling variability in NB models.
- **Null hypothesis (\(H_0\))**: No difference between groups for a gene.
- **p-value**: Probability of observing data as extreme as yours under \(H_0\).
- **Multiple testing**: Testing thousands of genes inflates false positives.
- **FDR (False Discovery Rate)**: Expected fraction of false positives among discoveries.
- **Benjamini–Hochberg (BH)**: Common procedure to control FDR.

---

## Differential Expression Outputs

- **log2 fold change (log2FC)**: Effect size; log2 ratio of expression between conditions.
- **Shrinkage**: Stabilizes noisy log2FC estimates (especially for low counts).
- **MA plot**: log2FC (M) vs mean expression (A) visualization.
- **Volcano plot**: log2FC vs \(-\log_{10}(p)\); highlights big effects with significance.
- **DEG**: Differentially expressed gene.

---

## Experimental Design Terms

- **Biological replicate**: Independent biological sample (e.g., different donors/animals).
- **Technical replicate**: Same sample processed multiple times; not a substitute for biology.
- **Confounder**: Variable correlated with condition (e.g., batch = condition).
- **Design matrix**: Encodes predictors (condition, batch, donor) for the model.
- **Contrast**: Specific comparison you test (e.g., treated vs control).

---

## Single-Cell Specific

- **Dropout / sparsity**: Many zero counts due to sampling and biology.
- **Pseudoreplication**: Treating cells as independent replicates when the true replicate is the sample/donor.
- **Pseudobulk**: Aggregating counts across cells per sample×celltype, then using bulk DE methods.

---

## Where You’ll Use These

- **Poisson/NB/dispersion**: Lab 2 + Module 8
- **p-values/FDR/BH**: Module 6 + Assignments 2–4
- **Design matrix/contrast/confounding**: Module 3 + Assignment 1 + DESeq2 labs
- **Pseudobulk/pseudoreplication**: Module 11 + Assignment 4


