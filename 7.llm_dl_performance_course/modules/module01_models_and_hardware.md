# Module 1: Models & Hardware

**Part 1 — Core inference**

---

## Learning objectives

By the end of this module you will be able to:

- Describe how LLMs and deep learning models interact with the hardware they run on.
- Distinguish inference from training in terms of compute and memory.
- Identify why “core pieces of models and hardware” are the right starting point for AI performance engineering.

---

## Summary

Efficient inference is one of the largest barriers to LLMs in production. To improve it, you need to understand both the **core pieces of the models** (e.g. layers, attention, tokenization, memory layout) and the **hardware** (GPUs, memory bandwidth, compute throughput). This module sets that foundation: how models map onto GPUs, where bottlenecks appear, and why inference differs from training (e.g. no gradients, different memory and batching patterns).

---

## Key concepts

- Inference vs training (compute and memory)
- Model structure (layers, parameters, activations) and memory
- GPU compute and memory bandwidth
- Why “study the core, then learn on demand” works for perf

---

## Curated resources (use these)

- **CS336 (Stanford), Lectures 5–11** — [link](https://buff.ly/9o4jRxl): model and systems aspects of LLMs.
- **Inside NVIDIA GPUs (A. Gordic)** — [link](https://buff.ly/wmgNwMg): how GPUs work under the hood.

---

## Connection to biology (Part 2)

The same ideas apply to gene-language models and cell-embedding models: their “core pieces” (layers, attention, embeddings) and the hardware (single GPU vs multi-GPU, memory limits) determine how fast you can run inference on large cell×gene datasets or serve models in production.
