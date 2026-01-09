# Glossary: Clustering & Cell Type Annotation

---

## Dimensionality Reduction

- **HVG (Highly Variable Genes)**: Genes with high variability across cells; used to focus on informative signal.
- **PCA (Principal Component Analysis)**: Linear projection capturing maximal variance directions.
- **PC loading**: Contribution of each gene to a principal component.
- **UMAP**: Non-linear embedding preserving local neighborhood structure (and some global).
- **t-SNE**: Non-linear embedding preserving local similarities; sensitive to perplexity/initialization.
- **n_neighbors (UMAP)**: Controls local vs global structure in UMAP.
- **perplexity (t-SNE)**: Effective number of neighbors; affects cluster separation/fragmentation.

---

## Clustering

- **k-NN graph**: Graph connecting each cell to its k nearest neighbors (in PCA/latent space).
- **SNN graph**: Graph weighting edges by shared neighbors; often stabilizes clustering.
- **Louvain / Leiden**: Community detection algorithms for graph-based clustering (Leiden improves connectivity).
- **Resolution**: Controls clustering granularity (higher → more clusters).
- **Stability**: Whether clusters persist across parameter choices / subsampling.
- **Silhouette score**: Measures how well-separated clusters are in embedding space (interpret with caution).

---

## Markers & Annotation

- **Marker gene**: Gene enriched in a cluster/cell type (ideally specific, not just high expression).
- **Specificity**: How unique a marker is to one cluster vs many.
- **Dotplot**: Shows expression fraction and average expression across clusters.
- **Reference atlas**: Curated annotated dataset used for label transfer.
- **SingleR**: Reference-based annotation method (R/Bioconductor).
- **CellTypist**: Marker/model-based annotation with confidence scoring (Python).
- **Confidence score**: Numeric score indicating how reliable an automated label is.
- **Granularity**: Level of detail in labels (e.g., “T cell” vs “CD4 memory T”).

---

## Artifacts

- **Doublet**: Two cells captured together; can masquerade as “new” cluster.
- **Technical covariate**: Non-biological driver (UMIs, mt%, ribosomal fraction).

---

## Where You’ll Use These

- **HVG/PCA/PC loadings**: `labs/lab02_pca.ipynb`
- **UMAP/t-SNE params**: `labs/lab03_umap_tsne.ipynb`
- **k-NN/Leiden/resolution**: `labs/lab05_clustering.ipynb`
- **Markers/dotplots**: `labs/lab07_markers.ipynb`
- **Automated annotation/confidence**: `labs/lab09_annotation.ipynb`


