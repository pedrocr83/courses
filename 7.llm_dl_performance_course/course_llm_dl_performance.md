# Course 7: LLMs, Deep Learning & Performance

## Efficient Inference and AI Performance Engineering at the Intersection with Biology

---

## Course Overview

| | |
|---|---|
| **Duration** | 3–4 weeks (15–22 hours total) |
| **Format** | Self-paced with curated resources and modules |
| **Level** | Intermediate (after Course 1; Courses 2/3 recommended) |
| **Prerequisites** | Course 1 (single-cell processing); Python, command line; basic idea of training vs inference |
| **Theme** | Core inference + biology application: AI/LLMs meet single-cell and omics |

## Start Here (Do This First)

- **Step-by-step checklist:** `START_HERE.md`
- **Glossary:** `glossary.md`
- **Environment setup:** `setup.md`
- **Curated readings:** `resources/reading/curated_readings.md`

### Learning Objectives

**Part 1 — Core inference**  
By the end of Part 1 you will be able to:

1. Explain why efficient inference is one of the largest barriers to LLMs in production.
2. Describe and reason about: batching, memory management, concurrent user sessions, GPU saturation, sharding and/or MIG (Multi-Instance GPUs), and parallelization.
3. Relate core model components and hardware (GPUs, memory bandwidth) to inference bottlenecks.
4. Use the curated courses, videos, and blogs as your primary references for going deeper.

**Part 2 — Biology application track**  
By the end of Part 2 you will be able to:

5. Apply Part 1 concepts to single-cell and omics workloads (e.g. gene-language models, cell embeddings, annotation at scale).
6. Reason about batching and memory for large cell×gene runs and model serving.
7. Identify where GPU saturation, MIG, or sharding matter in biological inference pipelines.

---

## Course Structure (Approach C)

### Part 1: Core LLM/DL & Inference

Study the core pieces of the models and the hardware they run on; then learn the rest on demand.

| Module | Topic | Key ideas |
|--------|--------|-----------|
| 1 | Models & hardware | How models and GPUs interact; inference vs training |
| 2 | Batching & memory | Batching, managing memory, impact on production |
| 3 | Concurrency & GPU saturation | Concurrent sessions, GPU utilization |
| 4 | Sharding & MIG | Sharding and Multi-Instance GPUs |
| 5 | Parallelization | Parallelization strategies for inference |

### Part 2: Biology Application Track

| Module | Topic | Key ideas |
|--------|--------|-----------|
| 6 | Single-cell & omics DL | Gene/cell models (e.g. scGPT, Geneformer); where inference and perf matter |
| 7 | Applying inference in biology | Batching, memory, scaling, serving; MIG/sharding for large atlases |

**On demand:** Use `resources/reading/curated_readings.md` for deeper dives; no fixed order.

---

## Curated Resources

All recommended courses, videos, blogs, and repos are listed in **`resources/reading/curated_readings.md`**, with short notes and mapping to Part 1 / Part 2. Use them as your main path to mastering AI performance engineering.

---

## Connection to the Rest of the Curriculum

- **Builds on:** Course 1 (processing, count matrices); Courses 2 and 3 help you connect inference to DE and clustering workflows.
- **Applies to:** Any later use of deep learning or LLMs in single-cell/omics (e.g. embedding models, gene-language models, annotation at scale).
