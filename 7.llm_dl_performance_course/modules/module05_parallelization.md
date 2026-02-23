# Module 5: Parallelization

**Part 1 — Core inference**

---

## Learning objectives

By the end of this module you will be able to:

- Describe main parallelization strategies for inference (data, tensor, pipeline).
- Explain how parallelization interacts with batching and memory.
- Reason about scaling inference across multiple GPUs or nodes.

---

## Summary

**Parallelization** means using multiple compute units (GPU cores, GPUs, or nodes) at once. For inference, common strategies include data parallelism (replicate model, split batch), tensor parallelism (split layers across devices), and pipeline parallelism (split layers in sequence). The right strategy depends on model size, batch size, and hardware.

---

## Key concepts

- Data vs tensor vs pipeline parallelism
- Communication and synchronization costs
- Scaling laws and bottlenecks

---

## Curated resources (use these)

- **AI Performance Engineering (C. Fregly)** — [link](https://buff.ly/8ey3GPa): performance and parallelization.
- **CS336 (Stanford), relevant lectures** — [link](https://buff.ly/9o4jRxl): systems and scaling.

---

## Connection to biology (Part 2)

When you scale gene or cell models to many GPUs (e.g. for a large atlas or high QPS), you’ll need to choose a parallelization strategy and possibly combine it with batching and MIG/sharding from earlier modules.
