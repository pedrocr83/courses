# Module 6: Single-Cell & Omics Deep Learning

**Part 2 — Biology application track**

---

## Learning objectives

By the end of this module you will be able to:

- Name representative gene-language and cell-embedding models (e.g. scGPT, Geneformer) and what they do.
- Identify where inference and performance matter in single-cell and omics pipelines.
- Connect Course 1–3 concepts (count matrices, clustering, annotation) to deep learning inference.

---

## Summary

Single-cell and omics increasingly use **deep learning**: gene-language models (trained on gene or transcript sequences), **cell-embedding models** (mapping cells to vectors), and **annotation or fate models**. These models are run at inference time on many cells or genes; batching, memory, and GPU use from Part 1 apply directly. This module grounds Part 1 in biology: what these models are and where their inference cost and throughput matter.

---

## Key concepts

- Gene-language models (e.g. Geneformer, scGPT) and cell-embedding models
- Inference in pipelines: embedding large atlases, annotation at scale
- Link to count matrices (Course 1), clustering (Course 3), and DE (Course 2)

---

## Curated resources (use these)

- Use **Part 1 resources** (NVIDIA, LLM Inference Patterns, etc.) and apply them to “model = gene/cell model.”
- Search for “scGPT,” “Geneformer,” “single-cell transformer” for biology-specific papers and repos.

---

## Connection to Part 1

Every concept from Part 1 applies: batching cells or sequences, memory for large models and long contexts, GPU saturation when embedding millions of cells, and (when needed) sharding or MIG for very large models or multi-user setups.
