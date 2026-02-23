# Module 7: Applying Inference in Biology

**Part 2 — Biology application track**

---

## Learning objectives

By the end of this module you will be able to:

- Apply batching and memory concepts to a single-cell or omics inference scenario.
- Reason about scaling and serving (e.g. embedding a large atlas, or serving an annotation API).
- Identify when MIG or sharding are relevant for biological inference workloads.

---

## Summary

This module ties Part 1 to practice: **batching** cells or sequences for throughput, **memory** limits and KV-cache–like behavior where relevant, **GPU saturation** when running embedding or annotation at scale, and **MIG/sharding** when the model or concurrency demands it. You’ll reason about “design a deployment” or “scale this pipeline” for a biology use case (e.g. cell-type annotation, gene-language model for variant effect, or atlas-wide embedding).

---

## Key concepts

- Batching and memory in a biology pipeline (cells, genes, requests)
- Scaling: single GPU → multi-GPU or MIG
- Serving: latency vs throughput for interactive vs batch use cases

---

## Curated resources (use these)

- **AI Performance Engineering (C. Fregly)** — [link](https://buff.ly/8ey3GPa): apply general perf concepts.
- **GPU Performance Engineering (E. Andere)** — [link](https://buff.ly/sfGJNMK): GPU-side tuning.
- Part 1 modules and `resources/reading/curated_readings.md` for reference.

---

## Connection to the curriculum

You’re now applying inference and performance (Course 7) to the same single-cell and omics workflows you learned in Courses 1–6: processing, DE, clustering, integration, trajectory, and cell-cell communication can all be combined with or augmented by DL models whose efficiency you can reason about.
