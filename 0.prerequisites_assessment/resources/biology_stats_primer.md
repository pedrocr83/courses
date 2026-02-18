# Biology & Statistics Primer for Single-Cell RNA-seq

**Estimated Time:** 10-15 hours (self-paced)
**Goal:** Build foundational knowledge to pass the Biology & Statistics Assessment (14/20)

---

## How to Use This Primer

- Work through both parts sequentially
- Take notes on key concepts—active recall improves retention
- Complete the review questions at the end before retaking the assessment
- Revisit sections where you feel uncertain

---

## Part 1: Molecular Biology Fundamentals

### 1. The Central Dogma: DNA → RNA → Protein

The **Central Dogma** describes the flow of genetic information in cells:

1. **DNA** stores the genetic blueprint. It is the permanent, heritable information.
2. **Transcription** copies a segment of DNA into RNA (DNA → RNA). This happens in the nucleus.
3. **Translation** converts the RNA message into a protein (RNA → Protein). This happens in the cytoplasm.

**DNA as the blueprint:** DNA is a double-stranded molecule. Each gene is a specific region of DNA that encodes instructions for making a product (usually a protein). DNA is not directly used to build proteins—it is first copied into RNA.

**Transcription (DNA → RNA):** An enzyme called RNA polymerase reads the DNA and synthesizes a complementary RNA strand. The RNA is a single-stranded copy of the gene. For protein-coding genes, this RNA is messenger RNA (mRNA).

**Translation (RNA → Protein):** Ribosomes read the mRNA and assemble amino acids into a polypeptide chain according to the genetic code. Each three-nucleotide codon specifies one amino acid. tRNA molecules ferry the correct amino acids to the ribosome. The resulting polypeptide folds into a functional protein.

**Why it matters for RNA-seq:** RNA-seq measures the RNA produced by transcription. We don't directly measure protein, but mRNA abundance is a strong predictor of protein levels for many genes (with caveats: mRNA stability, translation efficiency, and protein half-life all matter).

**Simple diagram (text):**
```
    [DNA]  ----transcription---->  [mRNA]  ----translation---->  [Protein]
  (nucleus)                       (nucleus                        (cytoplasm/
                                   to cytoplasm)                   ribosome)
```

**Key point:** Information flows from DNA → RNA → Protein. This unidirectional flow ensures that changes to RNA (e.g., from environmental stress) do not permanently alter the genome, while the genome remains the stable repository of genetic information.

---

### 2. What Is RNA?

**RNA** (ribonucleic acid) is a molecule similar to DNA but single-stranded and typically shorter-lived. Different types serve different functions:

| Type | Full Name | Function |
|------|-----------|----------|
| **mRNA** | Messenger RNA | Carries the code from DNA to ribosomes for protein synthesis |
| **rRNA** | Ribosomal RNA | Structural and catalytic component of ribosomes |
| **tRNA** | Transfer RNA | Brings amino acids to the ribosome during translation |
| **lncRNA** | Long non-coding RNA | Regulatory roles; does not code for proteins |
| **miRNA** | Micro RNA | Short RNAs that regulate gene expression by targeting mRNAs |

**mRNA as the messenger:** mRNA is the intermediary. A cell transcribes a gene → produces mRNA → translates mRNA into protein. The abundance of an mRNA reflects how actively that gene is being expressed.

**Why we care about mRNA in RNA-seq:** RNA-seq measures mRNA levels. By sequencing mRNA from cells, we infer which genes are expressed and at what levels. This tells us the cell's functional state. We typically capture mRNA (and sometimes other RNA types) and convert it to DNA for sequencing. mRNA makes up only a small fraction of total cellular RNA (rRNA dominates), so protocols often use poly(A) selection or rRNA depletion to enrich for mRNA.

---

### 3. Gene Expression

**Gene expression** means a gene is being used to produce its product (typically a protein). An expressed gene has been transcribed into RNA; if it codes for protein, that RNA may also be translated.

**Regulation:** Not all genes are active in all cells. Expression is controlled by transcription factors, epigenetic marks, and environmental signals. A liver cell expresses different genes than a neuron because each cell type has different functions.

**Expression levels:** Genes can be expressed at different levels:
- **High expression:** Many mRNA copies; the gene product is abundant.
- **Low expression:** Few mRNA copies; the gene product is rare.
- **Not expressed:** No detectable mRNA (gene is "off" or below detection).

Expression levels matter for function and for RNA-seq analysis—we must distinguish true low expression from technical noise or absence of expression.

**Upregulation and downregulation:** When we compare two conditions (e.g., treated vs control), a gene may be **upregulated** (higher expression) or **downregulated** (lower expression). These changes reflect altered regulation and often point to biological mechanisms.

---

### 4. Cell Types and Differentiation

**Same genome, different expression profiles:** Almost every cell in your body has the same DNA. What makes a neuron different from a muscle cell is *which genes are expressed* and *at what levels*—the gene expression profile.

**Stem cells → specialized cells:** Stem cells can divide and differentiate into specialized cell types. During differentiation, specific genes are turned on or off, leading to distinct expression programs and cell identities.

**Why different cell types exist:** Different functions require different proteins. A red blood cell needs hemoglobin; a pancreatic beta cell needs insulin. Each cell type has a characteristic set of expressed genes.

**Tissue heterogeneity:** Tissues are mixes of multiple cell types. A piece of liver contains hepatocytes, immune cells, endothelial cells, fibroblasts, and others. Single-cell RNA-seq reveals this heterogeneity by profiling each cell separately.

**Marker genes:** Genes that are specifically or preferentially expressed in one cell type are called marker genes. For example, CD3 is a T cell marker; hemoglobin genes mark erythrocytes. Marker genes help identify and annotate cell types in scRNA-seq.

---

### 5. The Transcriptome

**Genome vs transcriptome vs proteome:**
- **Genome:** The complete set of DNA in an organism. Same in (almost) all cells.
- **Transcriptome:** The complete set of RNA transcripts in a cell or population at a given time. Varies by cell type, condition, and time.
- **Proteome:** The complete set of proteins. Derived from the transcriptome but also influenced by translation efficiency and protein turnover.

**Why transcriptomics matters:** The transcriptome is a snapshot of which genes are active. It reflects developmental state, response to stimuli, and cell identity. It is easier to measure at scale than the proteome and more dynamic than the genome.

**Snapshot of cellular state:** RNA-seq captures the transcriptome at a moment in time. Cells change; repeating experiments at different times or under different conditions reveals how expression changes.

**Why not just measure protein?** Proteomics is advancing but remains more expensive and less sensitive for many applications. mRNA is an accessible proxy: it amplifies the signal (one gene can produce many mRNA copies) and is easier to capture at scale. For single-cell work, RNA-seq is currently the dominant technology.

---

### 6. Housekeeping Genes

**Housekeeping genes** are constitutively expressed in most or all cell types. They encode proteins needed for basic cellular functions (e.g., metabolism, structure, basic biosynthesis).

**Examples:**
- **GAPDH** (glyceraldehyde-3-phosphate dehydrogenase): Involved in glycolysis
- **ACTB** (β-actin): Structural protein, part of the cytoskeleton

**Why they matter as controls:** Because they are relatively stable across conditions and cell types, they are often used as:
- **Reference genes** for normalizing gene expression (e.g., qPCR, some RNA-seq methods)
- **Loading controls** in Western blots
- **Quality checks** in RNA-seq (low or absent housekeeping genes may indicate degradation)

Note: "Stable" is relative—housekeeping genes can vary in some contexts, so validation is important. In RNA-seq, global normalization methods (e.g., library size, DESeq2's median-of-ratios) have largely replaced single-gene references, but housekeeping genes remain useful for quality control.

---

### 7. Introduction to Single-Cell RNA-seq

**Bulk vs single-cell:** In **bulk RNA-seq**, RNA from many cells is pooled and sequenced together. You get average expression across the population. In **single-cell RNA-seq (scRNA-seq)**, each cell is profiled separately, so you see expression in individual cells.

**Why single-cell resolution matters:**
- Reveals cell type composition and rare cell types
- Captures heterogeneity within a population
- Allows tracking of developmental trajectories
- Identifies cell states and transitions

**Basic workflow overview (simplified):**
1. **Isolate single cells (or nuclei)** — dissociation, FACS, or microfluidics
2. **Capture and lyse each cell** — reverse-transcribe RNA to cDNA
3. **Amplify and add barcodes** — each molecule gets a cell barcode so we know which cell it came from
4. **Sequence the library** — high-throughput sequencing
5. **Align and quantify** — map reads to the reference genome, count per gene per cell

**What questions scRNA-seq can answer:**
- What cell types are in this tissue?
- How do cells transition from one state to another?
- Which genes define a cell type or state?
- How does gene expression vary across cells?
- What pathways or programs are active in each cell?

**Limitations to keep in mind:** Single-cell protocols capture only a fraction of each cell's mRNA (dropout), making some genes appear zero or low even when expressed. Sensitivity varies by protocol and cell size. These technical factors affect how we interpret the data.

---

## Part 2: Statistics Fundamentals

### 8. Descriptive Statistics

**Mean (average):** Sum of values divided by the number of values. Sensitive to outliers.

**Median:** The middle value when data are sorted. More robust to outliers than the mean.

**Mode:** The most frequently occurring value. Useful for categorical or discrete data.

**When to use each:**
- **Mean:** Symmetric data, no extreme outliers (e.g., many biological measurements).
- **Median:** Skewed data or when outliers are present.
- **Mode:** Categorical data (e.g., most common cell type).

**Standard deviation and variance:**
- **Variance** (σ²): Average squared deviation from the mean. Units are squared.
- **Standard deviation** (σ): Square root of variance. Same units as the data.

**Interpreting spread:** A large standard deviation means values are spread out; a small one means they are clustered near the mean. In RNA-seq, genes with high variance across cells are often biologically interesting (e.g., markers of cell states).

**Example:** If Gene A has mean expression 10 and std dev 2, and Gene B has mean 10 and std dev 8, Gene B varies much more across cells. That variation might reflect cell-type-specific expression or dynamic regulation.

---

### 9. Distributions

**Normal (Gaussian) distribution:** Bell-shaped, symmetric around the mean. Many biological and technical measurements approximate a normal distribution.

**Properties—68-95-99.7 rule:** For a normal distribution:
- ~68% of values fall within 1 standard deviation of the mean
- ~95% within 2 standard deviations
- ~99.7% within 3 standard deviations

**Skewed distributions:** Not symmetric. **Right-skewed:** long tail to the right (e.g., gene expression in RNA-seq, where most genes are low and few are very high). **Left-skewed:** long tail to the left.

**Histograms and how to read them:** Histograms show the distribution of a variable by binning values and showing counts or frequencies. Look at shape (symmetric vs skewed), center, spread, and outliers.

**Why RNA-seq data are often log-transformed:** Gene expression counts are typically right-skewed (many low values, few very high). Log transformation (e.g., log2(counts+1)) makes distributions more symmetric and stabilizes variance, which improves the behavior of many statistical methods and visualizations.

---

### 10. Hypothesis Testing

**Null hypothesis (H₀):** The default assumption (e.g., "no difference," "no effect"). We test whether the data provide evidence against it.

**Alternative hypothesis (H₁):** The opposite of the null (e.g., "there is a difference").

**p-values:** The probability of observing data as extreme as (or more extreme than) what we saw, *if the null hypothesis were true*. A small p-value suggests the null is unlikely.

**What p-values do NOT mean:**
- They are NOT the probability that the null is true
- They are NOT the probability that a result is due to chance
- They do NOT measure effect size or biological importance

**Significance threshold (α = 0.05):** Conventionally, we reject the null if p < 0.05. This means we accept up to a 5% chance of a false positive when the null is true.

**Type I errors (false positives):** Rejecting the null when it is actually true. Probability = α.

**Type II errors (false negatives):** Failing to reject the null when it is actually false. Related to statistical power (ability to detect true effects). Low power (e.g., too few replicates) increases the risk of missing real differences.

**Example in RNA-seq:** We test "Is Gene X differentially expressed between condition A and B?" H₀: no difference. If p = 0.02, we reject H₀ and conclude there is evidence of a difference. If p = 0.15, we do not reject H₀—we lack evidence, but that doesn't prove the gene is unchanged.

---

### 11. Multiple Testing Problem

**The problem:** If you test 20,000 genes at α = 0.05, you expect 20,000 × 0.05 = 1,000 false positives *even when no gene is truly different*. With many tests, unadjusted p-values are misleading.

**Family-wise error rate (FWER):** Probability of making at least one Type I error across all tests. Bonferroni correction controls FWER by using α/n as the per-test threshold (very strict for large n).

**False Discovery Rate (FDR):** The expected proportion of significant results that are false positives. More lenient than FWER when many true positives are expected.

**Benjamini-Hochberg correction:** A common method to control FDR. It adjusts p-values so that when you call results significant at FDR < 0.05, about 5% of those are expected to be false discoveries. The procedure: sort p-values, find the largest rank *k* where pₖ ≤ (k/n)×α, and reject hypotheses 1 through *k*.

**Practical takeaway:** In RNA-seq, we almost always use FDR or similar multiple-testing corrections. A raw p-value of 0.01 among 20,000 genes is not compelling without correction. Software like DESeq2, edgeR, and Scanpy perform these corrections automatically—look for "adjusted p-value," "padj," or "q-value" in outputs.

---

### 12. Correlation

**Pearson correlation (r):** Measures linear association between two continuous variables. Ranges from -1 to 1. Assumes normality.

**Spearman correlation (ρ):** Rank-based; measures monotonic association. More robust to outliers and non-normal distributions.

**Correlation does NOT imply causation:** Two variables can be correlated because one causes the other, because a third variable causes both, or by chance. Correlation alone cannot establish causation.

**Examples from biology:** Gene A and Gene B may be correlated because they are co-regulated, because one regulates the other, or because they respond to the same stimulus.

**Interpreting r values:**
- |r| close to 1: Strong linear relationship
- |r| around 0.5: Moderate relationship
- |r| close to 0: Weak or no linear relationship

**In single-cell analysis:** We often compute correlation between genes across cells (e.g., to find co-expressed genes or build gene modules) or between cells (e.g., for clustering). Spearman is commonly used because expression data are often not normally distributed.

---

### 13. Statistical Significance vs Biological Significance

**A gene can be statistically significant but biologically irrelevant:** With enough samples or very precise measurements, tiny differences can have very small p-values. A 0.1% change in expression might be statistically significant but biologically meaningless.

**Effect size matters:** In RNA-seq, we care about **log2 fold change**. A gene with log2FC = 2 (4-fold change) is typically more biologically interesting than one with log2FC = 0.1. Always consider both p-value and effect size.

**Practical significance in genomics:** We look for genes that are both statistically significant (after multiple-testing correction) and have meaningful effect sizes (e.g., |log2FC| > 1 for "large" effects). Biological context and prior knowledge also matter.

**Log2 fold change recap:** log2FC = 1 means 2-fold change; log2FC = 2 means 4-fold; log2FC = -1 means half the expression. It is used because fold changes span many orders of magnitude, and log scale makes comparisons symmetric (e.g., 2× up and 0.5× down are equidistant from 0).

---

## Review Questions

### Biology (5 questions)

**Q1.** What are the two main steps of the Central Dogma that produce protein from DNA, and where does each occur?

**A1.** Transcription (DNA → RNA) occurs in the nucleus. Translation (RNA → Protein) occurs in the cytoplasm at ribosomes.

---

**Q2.** Why do we focus on mRNA in RNA-seq instead of other RNA types?

**A2.** mRNA carries the coding information from genes to ribosomes. Its abundance reflects how actively each gene is being expressed, giving a snapshot of the cell's functional state. rRNA and tRNA are more abundant but less informative for which genes are active.

---

**Q3.** Why can a liver cell and a neuron have different functions if they have the same genome?

**A3.** Both have the same DNA, but they express different sets of genes at different levels. The gene expression profile defines cell identity and function. Different transcription factors and epigenetic states lead to different expression programs.

---

**Q4.** What is the difference between the genome and the transcriptome?

**A4.** The genome is the complete set of DNA (largely static and shared across cells). The transcriptome is the complete set of RNA molecules present at a given time (dynamic, varies by cell type, condition, and time).

---

**Q5.** What advantage does single-cell RNA-seq have over bulk RNA-seq?

**A5.** Single-cell RNA-seq profiles each cell individually, revealing cell type composition, heterogeneity within populations, rare cell types, and cell-state transitions. Bulk RNA-seq only gives average expression across a mixed population.

---

### Statistics (5 questions)

**Q6.** When would you use the median instead of the mean to describe the center of a distribution?

**A6.** Use the median when the distribution is skewed or when there are influential outliers. The median is more robust because it is not pulled by extreme values.

---

**Q7.** What does a p-value of 0.03 mean? What does it NOT mean?

**A7.** It means that if the null hypothesis were true, we would see data as extreme as ours (or more extreme) 3% of the time. It does NOT mean there is a 3% chance the null is true, nor that the result has a 97% chance of being correct.

---

**Q8.** Why is it important to correct for multiple testing when analyzing thousands of genes?

**A8.** With thousands of tests at α = 0.05, we expect many false positives by chance. Correction (e.g., Benjamini-Hochberg for FDR) controls the rate of false discoveries and makes our list of significant genes more reliable.

---

**Q9.** If two genes have a Pearson correlation of 0.9, can we conclude that one causes the other? Why or why not?

**A9.** No. Correlation measures association, not causation. The genes could be co-regulated, one could regulate the other, or both could be driven by a third factor. Correlation alone cannot establish causality.

---

**Q10.** A gene has an adjusted p-value of 0.001 but a log2 fold change of 0.05. Should it be considered biologically important? Why or why not?

**A10.** Probably not. It is statistically significant but the effect size is very small (≈3.5% change). Such a difference is often within technical noise and rarely has meaningful biological impact. Both statistical significance and effect size should be considered.

---

## Additional Resources

### Molecular Biology & Genomics
- [Khan Academy: Gene Expression (Central Dogma)](https://www.khanacademy.org/science/biology/gene-expression-central-dogma) — Free video lessons
- [NHGRI: Introduction to Genomics](https://www.genome.gov/About-Genomics/Introduction-to-Genomics) — NIH overview
- [DNA Learning Center](https://dnalc.cshl.edu/) — Interactive modules

### Statistics
- [Khan Academy: Statistics & Probability](https://www.khanacademy.org/math/statistics-probability) — Comprehensive free course
- [Seeing Theory](https://seeing-theory.brown.edu/) — Visual, interactive statistics
- [StatQuest (YouTube)](https://www.youtube.com/user/joshstarmer) — Clear explanations of statistical concepts

### RNA-seq and Single-Cell
- [Harvard Bioinformatics Core (HBC) Training](https://hbctraining.github.io/) — RNA-seq tutorials
- [Single-cell best practices](https://www.sc-best-practices.org/) — Community-driven scRNA-seq workflows

### Supplemental
- [StatQuest: RNA-seq](https://www.youtube.com/playlist?list=PLblh5JKOoLUJ5fjyQp6fDExKQJgLfk6hR) — Accessible videos on RNA-seq concepts
- [Genomics 101 (NIH)](https://www.genome.gov/genetics-glossary) — Glossary of genomic terms

---

## Study Tips

1. **Spaced repetition:** Review sections 2–3 days after first reading.
2. **Draw it:** Sketch the Central Dogma, a cell differentiation tree, or a normal distribution.
3. **Connect to single-cell:** As you learn each concept, ask how it applies to scRNA-seq (e.g., "Why would we log-transform counts?").
4. **Practice the math:** Compute mean, median, and std dev for a small dataset; interpret a correlation matrix.
5. **Test yourself:** Cover the answers to the review questions and try answering before peeking—active recall strengthens memory.

---

## Part 1 Summary: Biology Concepts for scRNA-seq

| Concept | Why It Matters |
|---------|----------------|
| Central Dogma | RNA-seq measures the transcript; understand the pathway from gene to protein |
| mRNA | The primary target of RNA-seq; reflects gene expression |
| Gene expression & regulation | Cells differ by which genes are on/off and at what levels |
| Cell types & markers | scRNA-seq identifies cell types via expression profiles and marker genes |
| Transcriptome | RNA-seq captures this; it's dynamic and cell-type-specific |
| Housekeeping genes | Used for QC and (historically) normalization |
| Bulk vs single-cell | Single-cell reveals heterogeneity; bulk averages over it |

---

## Part 2 Summary: Statistics Concepts for scRNA-seq

| Concept | Why It Matters |
|---------|----------------|
| Mean, median, std dev | Describe expression distributions; median robust to outliers |
| Distributions & log transform | Expression data are skewed; log stabilizes variance |
| Hypothesis testing & p-values | We test whether genes differ between conditions |
| Multiple testing & FDR | Correct for testing many genes; use adjusted p-values |
| Correlation | Find co-expressed genes; Spearman for non-normal data |
| Effect size (log2FC) | Distinguish biologically meaningful from trivial differences |

---

## Before You Retake the Assessment

- Review any sections where you hesitated on the practice questions.
- Ensure you can explain the Central Dogma, gene expression, and the transcriptome in your own words.
- Be comfortable interpreting p-values, FDR, and log2 fold change.
- The assessment has 20 questions (biology and statistics); passing requires 14 correct (70%).

*After completing this primer and the review questions, retake the Biology & Statistics Assessment. Aim for 14/20 or higher to proceed.*
