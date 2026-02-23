# Module 2: Batching & Memory

**Part 1 — Core inference**

---

## Learning objectives

By the end of this module you will be able to:

- Explain what batching is and why it matters for inference throughput and GPU utilization.
- Describe the main memory considerations in inference (weights, activations, KV cache).
- Relate batching and memory management to production deployment.

---

## Summary

**Batching** means processing multiple inputs (e.g. requests, sequences) together so the GPU does more useful work per kernel launch and achieves better saturation. **Memory management** involves fitting model weights, activations, and (in autoregressive LLMs) KV caches in device memory, and moving data efficiently. Together, batching and memory are central to making inference efficient and predictable in production.

---

## Key concepts

- Batching: dynamic vs static batching; batch size vs latency
- Memory: weights, activations, KV cache; memory bandwidth
- Trade-offs: throughput vs latency, batch size vs OOM

---

## Curated resources (use these)

- **LLM Inference Patterns** — [link](https://buff.ly/4ASDS5c): patterns that include batching and memory.
- **Understanding LLM Inference (article)** — [link](https://buff.ly/NrNU3Qz): deep dive on how inference works and where memory matters.

---

## Connection to biology (Part 2)

In single-cell and omics, “inputs” might be cells or genes: batching cells for embedding, or batching gene sequences for a gene-language model. Memory limits how many cells or how long a sequence you can process at once; batching determines throughput when you scale to large atlases or high request rates.
