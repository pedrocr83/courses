# Glossary: Single-Cell RNA-seq Sample Processing

Use this as a quick reference while doing labs and assignments.

---

## Core Terms

- **FASTQ**: Text format for sequencing reads; each read has 4 lines (header, sequence, separator, quality).
- **Read 1 / Read 2 (R1/R2)**: Paired-end reads; in 10x, R1 often contains barcode/UMI and R2 contains cDNA.
- **Cell barcode**: Short sequence identifying the droplet/cell of origin.
- **UMI (Unique Molecular Identifier)**: Short sequence used to deduplicate PCR duplicates and estimate unique molecules.
- **cDNA**: Complementary DNA synthesized from RNA; sequenced to infer transcripts.
- **Phred quality score**: Encodes base call accuracy; higher = lower error probability.
- **Reference genome / transcriptome**: FASTA + GTF annotation used for alignment/quantification.
- **GTF**: Gene annotation file specifying exon/gene/transcript coordinates.
- **Alignment**: Mapping reads to reference (e.g., STAR/Cell Ranger).
- **Pseudo-alignment**: Faster mapping to transcripts without full alignment (e.g., kallisto).
- **Count matrix**: Genes × cells matrix of UMI counts.
- **Sparse matrix**: Efficient representation when most entries are zeros (common in scRNA-seq).

---

## Cell Calling & Ambient RNA

- **Cell calling**: Distinguishing true cells from empty droplets based on counts/barcode rank.
- **Barcode rank (knee) plot**: Counts per barcode vs rank on log-log; knee/inflection suggests cell boundary.
- **Empty droplets**: Droplets without cells that still capture ambient RNA.
- **Ambient RNA**: Free-floating RNA captured by droplets, contaminating counts.

---

## Quality Control Metrics

- **`total_counts`**: Total UMIs per cell (library size).
- **`n_genes_by_counts`**: Number of genes detected per cell.
- **`pct_counts_mt`**: Percent of counts from mitochondrial genes; often higher in dying cells.
- **Outlier filtering**: Removing cells with metrics far from typical distribution (thresholds or MAD-based).
- **MAD (Median Absolute Deviation)**: Robust dispersion measure for outlier detection.

---

## Artifacts

- **Doublet**: Two cells captured in one droplet; appears like hybrid expression profile.
- **Scrublet / DoubletFinder**: Tools to score/predict doublets.
- **PCR duplicate**: Same original molecule amplified; UMIs help deduplicate.

---

## Outputs & Reproducibility

- **H5AD**: AnnData storage format used by Scanpy.
- **MTX**: Matrix Market format for sparse matrices.
- **Reproducible pipeline**: Re-runnable workflow with recorded versions/parameters and no manual “mystery steps”.

---

## Where You’ll Use These

- **FASTQ, barcode, UMI, Phred**: Lab 3 + Assignment 1
- **Knee plot, Empty droplets, Ambient RNA**: Lab 6
- **QC metrics, MAD, filtering**: Lab 8–9 + Assignment 3
- **Doublets**: Lab 9


