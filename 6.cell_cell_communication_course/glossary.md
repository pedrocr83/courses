# Glossary: Cell–Cell Communication (CCC)

---

## Biology

- **Ligand**: Secreted or membrane-bound molecule that binds a receptor to trigger signaling.
- **Receptor**: Protein (often membrane) that binds ligand; initiates downstream signaling.
- **Co-receptor**: Additional receptor component required for signaling.
- **Autocrine / paracrine / juxtacrine / endocrine**: Signaling to self / nearby / contact-dependent / long-range.
- **Pathway**: Set of interactions in a signaling cascade (e.g., TGFβ, Notch).

---

## Inference Concepts

- **L–R pair**: A ligand–receptor interaction from a curated database.
- **Interaction score**: A numeric score quantifying evidence for a given sender→receiver interaction.
- **Permutation test**: Statistical test using label shuffling (CellPhoneDB core approach).
- **Dropout**: Zeros due to sampling; can hide true expression.
- **Thresholding**: Minimum expression/percent-expressed rules to reduce noise.

---

## Tools

- **CellPhoneDB**: Python tool; permutation-based CCC inference using curated L–R pairs.
- **CellChat**: R tool; infers communication probabilities and summarizes at pathway level.
- **NicheNet**: R tool; links sender ligands to receiver target gene programs (ligand activity).
- **LIANA**: Framework (R/Python) to run multiple methods and compute consensus rankings.

---

## Networks

- **Graph / network**: Nodes = cell types; edges = interactions; weights = strength/score.
- **Centrality**: Network importance metric (sender/receiver “hubness”).
- **Community**: Subnetworks/modules of strongly interacting cell types.

---

## Validation

- **Orthogonal validation**: Independent evidence (protein, imaging, perturbation) supporting predicted interactions.
- **Spatial constraint**: Communication requires proximity for many signals; spatial data can confirm plausibility.

---

## Where You’ll Use These

- **Databases and coverage**: `labs/lab01_lr_databases.Rmd`
- **Permutation testing / CellPhoneDB outputs**: `labs/lab04_cellphonedb.Rmd`
- **Pathway-level inference / centrality**: `labs/lab05_cellchat.Rmd`
- **Consensus ranking**: `labs/lab07_liana.Rmd`
- **Network analysis**: `labs/lab10_network_analysis.Rmd`


