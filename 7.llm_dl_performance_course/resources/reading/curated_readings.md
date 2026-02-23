# Curated Readings: LLMs, Deep Learning & Performance (Course 7)

Use these as your main path to mastering AI performance engineering. Study the core pieces of the models and the hardware they run on; then learn the rest on demand.

---

## 1. Courses / Videos

| Topic | Title | Link | Part | Notes |
|-------|--------|------|------|--------|
| LLM foundations & systems | CS336 Stanford (Lectures 5–11) | [CS336](https://buff.ly/9o4jRxl) | Part 1 | Model and systems aspects of LLMs |
| GPU programming & perf | GPU MODE | [GPU MODE](https://buff.ly/bPEHYL0) | Part 1 | How to think about GPU utilization |
| Inference patterns | LLM Inference Patterns | [Inference Patterns](https://buff.ly/4ASDS5c) | Part 1 | Batching, memory, serving patterns |
| Production inference | LLM Inference (NVIDIA) | [NVIDIA LLM Inference](https://buff.ly/muDMcor) | Part 1 | Concurrency, MIG, scaling |

---

## 2. Repositories / Blogs

| Topic | Title | Link | Part | Notes |
|-------|--------|------|------|--------|
| GPU perf engineering | GPU Performance Engineering (E. Andere) | [Andere](https://buff.ly/sfGJNMK) | Part 1 | Sharding, MIG, tuning |
| AI perf engineering | AI Performance Engineering (C. Fregly) | [Fregly](https://buff.ly/8ey3GPa) | Part 1, Part 2 | Parallelization, scaling |
| Inside GPUs | Inside NVIDIA GPUs (A. Gordic) | [Gordic](https://buff.ly/wmgNwMg) | Part 1 | Hardware and architecture |
| Understanding inference | Understanding LLM Inference (article) | [Article](https://buff.ly/NrNU3Qz) | Part 1 | Batching, memory, how inference works |

---

## How to use this list

- **Part 1 (core inference):** Start with CS336 (relevant lectures), “Understanding LLM Inference,” and LLM Inference Patterns. Add GPU MODE and NVIDIA for hardware and production. Use Andere and Fregly for deep dives on GPU and parallelization.
- **Part 2 (biology):** Apply the same resources to gene/cell models and single-cell pipelines; no separate “biology-only” list—the intersection is applying Part 1 to biology.
- **On demand:** There is no single step-by-step guide to “master” AI perf engineering; the best bet is to study the core model and hardware pieces, then learn the rest on demand from these resources.
