# Assignment 2: Biology Application (Part 2)

**Complete after Part 2 (Modules 6–7).**

---

## Objective

Apply Part 1 inference concepts to a single-cell or omics scenario: batching, memory, scaling, and (if relevant) MIG or sharding.

---

## Task

**Design a small deployment scenario** for a biology inference workload. Choose one of (or combine):

- **Option A:** Serving a cell-embedding model (e.g. for cell-type annotation or atlas embedding) with multiple users or batch jobs.
- **Option B:** Running a gene-language model on many sequences (e.g. variant effect or gene representation) with limited GPU memory.

In 1–2 pages (or equivalent), address:

1. **Workload:** What model, what inputs (cells/genes/requests), what throughput or latency goal (rough).
2. **Batching:** How you would batch inputs and why.
3. **Memory:** What you’d watch (model size, batch size, OOM risk).
4. **Scaling:** Would you use a single GPU, MIG, or multi-GPU sharding? Why?
5. **One trade-off:** e.g. latency vs throughput, or cost vs simplicity.

---

## Deliverable

- A short design document (≈1–2 pages) saved to your `results/` or `notes/` folder (e.g. `assignment2_biology_application.md` or `.pdf`).

---

## Self-check

- [ ] Workload is clearly described (model + inputs + goal).
- [ ] At least batching and one of memory/scaling/MIG/sharding are discussed.
- [ ] Trade-off is stated in one or two sentences.
