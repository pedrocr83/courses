# Module 4: Sharding & MIG

**Part 1 — Core inference**

---

## Learning objectives

By the end of this module you will be able to:

- Explain what sharding is (splitting model or workload across devices) and when it is used.
- Describe MIG (Multi-Instance GPU): what it is, when to use it, and how it differs from full-GPU scheduling.
- Choose between “one big GPU,” MIG, and multi-GPU sharding for a given scenario.

---

## Summary

**Sharding** splits a model or workload across multiple GPUs (or nodes) to fit in memory or increase throughput. **MIG** (NVIDIA) partitions a single GPU into smaller instances so multiple jobs or users can share the GPU with isolation. Both are tools for scaling and multi-tenancy; the right choice depends on model size, request patterns, and cost.

---

## Key concepts

- Model parallelism vs data parallelism; tensor and pipeline parallelism
- MIG: partitions, isolation, use cases
- When to shard vs when to use MIG vs when one GPU is enough

---

## Curated resources (use these)

- **NVIDIA LLM Inference** — [link](https://buff.ly/muDMcor): MIG and multi-GPU inference.
- **GPU Performance Engineering (E. Andere)** — [link](https://buff.ly/sfGJNMK): GPU perf and partitioning.

---

## Connection to biology (Part 2)

Large gene-language models or cell-embedding models may not fit on one GPU; sharding lets you run them. MIG can let several researchers or jobs share a single GPU (e.g. different atlases or batch sizes) with predictable performance and isolation.
