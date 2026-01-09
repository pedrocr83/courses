# Glossary: Data Integration & Batch Correction

---

## Core Concepts

- **Batch effect**: Technical differences between datasets/samples that distort biological structure.
- **Confounding**: When batch and condition are correlated (e.g., all controls in batch 1, all treated in batch 2).
- **Integration**: Aligning datasets into a shared space while preserving biology.
- **Batch correction**: Removing unwanted technical variation (may be done in expression space or embedding space).

---

## Methods

- **ComBat**: Empirical Bayes batch correction method (often used for bulk/pseudobulk).
- **Regression correction**: Modeling and removing covariate effects (use carefully to avoid removing biology).
- **MNN (Mutual Nearest Neighbors)**: Uses pairs of nearest neighbors across batches to estimate correction vectors.
- **Harmony**: Corrects PCA embeddings via iterative soft clustering to mix batches.
- **Seurat CCA/RPCA**: Anchor-based integration; CCA for similar datasets, RPCA for mapping queries to references.
- **scVI**: Variational autoencoder producing a batch-corrected latent space.
- **scANVI**: Semi-supervised variant of scVI that uses labels to improve integration.

---

## Evaluation

- **Batch mixing**: Whether cells from different batches overlap within the same biology.
- **Biological conservation**: Whether true cell types/states remain distinct after correction.
- **LISI (Local Inverse Simpson Index)**: Diversity score in local neighborhoods.
  - **iLISI**: batch mixing (higher is better).
  - **cLISI**: cell-type conservation (lower is better).
- **kBET**: Tests whether local neighborhoods reflect expected batch proportions.
- **ARI (Adjusted Rand Index)**: Agreement between cluster assignments and known labels.
- **ASW (Average Silhouette Width)**: Separation measure; can be computed for batch or cell types.

---

## Practical

- **Reference mapping**: Project query dataset onto a well-annotated reference; often safer than full integration.
- **Over-correction**: Removing real biological differences.
- **Under-correction**: Leaving batch effects unresolved.

---

## Where You’ll Use These

- **Batch/confounding diagnosis**: `labs/lab02_batch_diagnosis.ipynb`
- **Harmony**: `labs/lab05_harmony.ipynb`
- **scVI**: `labs/lab07_scvi.ipynb`
- **Evaluation metrics**: `labs/lab08_evaluation.ipynb`
- **Reference mapping**: `labs/lab10_reference_mapping.ipynb`


