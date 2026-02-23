# Module 3: Concurrency & GPU Saturation

**Part 1 — Core inference**

---

## Learning objectives

By the end of this module you will be able to:

- Explain what concurrent user sessions mean for inference systems and how they relate to batching.
- Define GPU saturation and why it matters for cost and throughput.
- Describe how to measure and improve GPU utilization in practice.

---

## Summary

**Concurrent user sessions** mean multiple users or requests are being served at the same time. To do that efficiently, you need batching, scheduling, and often multiple GPUs or MIG. **GPU saturation** means keeping the GPU busy (high utilization of compute and memory bandwidth) instead of idle. This module links concurrency, batching, and saturation so you can reason about real-world serving and multi-tenant setups.

---

## Key concepts

- Concurrent requests and sessions; scheduling
- GPU utilization and saturation; why under-utilization is costly
- Measuring utilization (e.g. nvidia-smi, profilers)

---

## Curated resources (use these)

- **NVIDIA LLM Inference** — [link](https://buff.ly/muDMcor): production inference and utilization.
- **GPU MODE** — [link](https://buff.ly/bPEHYL0): GPU programming and performance.

---

## Connection to biology (Part 2)

In biology pipelines, “concurrent sessions” might be multiple analysts or jobs (e.g. embedding different atlases, or running cell-type annotation at scale). Saturation matters when you pay for GPU time: better batching and scheduling mean more cells or requests per hour per GPU.
