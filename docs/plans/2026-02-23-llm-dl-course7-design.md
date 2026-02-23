# Course 7: LLMs, Deep Learning & Performance — Design (Approach C)

**Approach:** Core inference + biology application track.  
**Theme:** Intersection of AI/LLMs with biology and single-cell (same program as Courses 1–6).

---

## 1. Goal & Scope

- **Goal:** Add Course 7 to the existing single-cell curriculum so learners master **efficient LLM/deep-learning inference** (batching, memory, concurrency, GPU saturation, sharding/MIG, parallelization) and apply it in **single-cell/biology** contexts.
- **Scope:** Two-part structure: **Part 1** = core models + hardware + inference (using curated courses/videos and blogs); **Part 2** = biology application track (single-cell/omics: gene/cell models, scaling, deployment).

---

## 2. Placement & Prerequisites

- **Placement:** Course 7 in the same program; listed in root README and docs like Courses 1–6.
- **Prerequisites:**
  - **From curriculum:** Course 1 (single-cell processing) completed; Course 2 or 3 recommended so learners know count matrices, DE, clustering.
  - **General:** Python, command line, basic idea of “training vs inference”; no GPU experience required to start.
- **Path:** “Core pieces of models and hardware first; then application to single-cell; rest on demand.”

---

## 3. Learning Objectives

**Part 1 (Core inference)**  
By the end of Part 1, learners will be able to:

- Explain why efficient inference is a major barrier to LLMs in production.
- Describe and use: batching, memory management, concurrent user sessions, GPU saturation, sharding and/or MIG, parallelization.
- Map core model components and hardware (GPUs, memory bandwidth) to inference bottlenecks.
- Use curated resources (CS336, GPU MODE, NVIDIA LLM Inference, inference-patterns, blogs) as primary references.

**Part 2 (Biology application track)**  
By the end of Part 2, learners will be able to:

- Apply Part 1 concepts to single-cell/omics workloads (e.g. gene-language models, cell embeddings, annotation at scale).
- Reason about batching and memory for large cell × gene runs and model serving.
- Identify where GPU saturation, MIG, or sharding matter in biological inference pipelines.

---

## 4. Module Outline (Approach C)

### Part 1: Core LLM/DL & Inference

| Module | Topic | Content (concise) | Curated resources |
|--------|--------|-------------------|-------------------|
| 1 | Models & hardware | Core pieces of models and hardware they run on; inference vs training | CS336 (e.g. 5–11), Inside NVIDIA GPUs |
| 2 | Batching & memory | Batching, managing memory, when they matter in production | LLM Inference Patterns, Understanding LLM Inference (article) |
| 3 | Concurrency & GPU saturation | Concurrent user sessions, GPU saturation, measuring utilization | NVIDIA LLM Inference, GPU MODE |
| 4 | Sharding & MIG | Sharding and/or Multi-Instance GPUs; when to use | NVIDIA, GPU Performance Engineering (E. Andere) |
| 5 | Parallelization | Parallelization strategies for inference | AI Performance Engineering (C. Fregly), CS336 |

### Part 2: Biology Application Track

| Module | Topic | Content (concise) |
|--------|--------|-------------------|
| 6 | Single-cell & omics DL | Gene/cell models (e.g. scGPT, Geneformer); where inference and perf matter |
| 7 | Applying inference concepts | Batching and memory for cell/gene models; scaling and serving in practice; optional: MIG/sharding for large atlases |

**On-demand:** Deeper dives via curated readings and links; no fixed order.

---

## 5. Course Structure (Repo)

- **Folder:** `7.llm_dl_performance_course/`
- **Same conventions as other courses:** `START_HERE.md`, `course_llm_dl_performance.md`, `setup.md`, `glossary.md`, `modules/`, `resources/` (including `resources/reading/curated_readings.md`), `quizzes/`, `assignments/`, `labs/` (optional/light in v1).
- **Curated readings:** User’s list (courses, videos, blogs, repos) in `resources/reading/curated_readings.md` with short notes and mapping to Part 1 / Part 2.

---

## 6. Integration with Existing Program

- **README:** Add Course 7 section (title, duration, location, START_HERE link, 1–2 sentence description, “Part 1: Core inference / Part 2: Biology application”).
- **docs/MASTER_TRAINING_GUIDE.md:** Add Course 7 to curriculum structure and prerequisite graph (after Course 1; Courses 2/3 recommended).
- **config/notebooklm_courses.yaml:** Add entry for `7.llm_dl_performance_course` with same pattern as other courses (START_HERE, course doc, modules, resources/reading, quizzes).
- **docs/ASSIGNMENT_RUBRICS.md:** Add placeholders for Course 7 assignments if rubrics are used.

---

## 7. Success Criteria

- Course 7 is discoverable and consistent with Courses 1–6 (structure, START_HERE, prerequisites).
- Part 1 is clearly “core inference”; Part 2 is clearly “biology/single-cell application.”
- All user-provided resources are linked and briefly described in `curated_readings.md`.
- A learner who completes Part 1 + Part 2 can explain inference bottlenecks and apply concepts to single-cell/omics scenarios.

---

**Design approved; implementation plan follows in `docs/plans/2026-02-23-llm-dl-course7-implementation.md`.**
