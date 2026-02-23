# Course 7: LLMs, Deep Learning & Performance — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add Course 7 (Approach C: core inference + biology application track) to the single-cell curriculum with full folder structure, modules, curated readings, and repo integration.

**Architecture:** New course directory `7.llm_dl_performance_course/` following existing course conventions; Part 1 modules = core inference; Part 2 modules = biology application; user resources in `resources/reading/curated_readings.md`; README and docs updated.

**Tech Stack:** Markdown, YAML (config); no new runtimes.

---

## Task 1: Create course directory and core files

**Files:**
- Create: `7.llm_dl_performance_course/`
- Create: `7.llm_dl_performance_course/course_llm_dl_performance.md`
- Create: `7.llm_dl_performance_course/setup.md`
- Create: `7.llm_dl_performance_course/glossary.md`

**Step 1:** Create directory `7.llm_dl_performance_course/`.

**Step 2:** Create `course_llm_dl_performance.md` with: Course Overview table (Duration, Format, Level, Prerequisites); Start Here / Glossary / Setup pointers; Learning Objectives for Part 1 and Part 2; Course Structure (Part 1 modules 1–5, Part 2 modules 6–7); references to `resources/reading/curated_readings.md`.

**Step 3:** Create `setup.md` with: Python 3.10+; optional GPU (CUDA) for labs; no strict R dependency; links to CS336, GPU MODE, NVIDIA if relevant for setup.

**Step 4:** Create `glossary.md` with: Batching, MIG, sharding, GPU saturation, inference (vs training), parallelization, and 5–10 other terms used in Part 1/2.

**Step 5:** Commit: `git add 7.llm_dl_performance_course/ && git commit -m "feat(course7): add course outline, setup, glossary"`

---

## Task 2: Add START_HERE.md

**Files:**
- Create: `7.llm_dl_performance_course/START_HERE.md`

**Step 1:** Copy structure from `docs/STANDARDIZED_START_HERE_TEMPLATE.md` and adapt for Course 7: orientation (course_llm_dl_performance.md, glossary, resources); setup (setup.md); learning calendar (Part 1 weeks → Part 2 weeks with module/lab refs); quizzes and assignments sections; link to curated readings; progress tracking; “Before moving to next course” (optional “next” or end of program).

**Step 2:** Ensure all paths point to files that will exist (e.g. `modules/module01_*.md`, `quizzes/quiz_module01.md`).

**Step 3:** Commit: `git add 7.llm_dl_performance_course/START_HERE.md && git commit -m "feat(course7): add START_HERE checklist"`

---

## Task 3: Create Part 1 modules (core inference)

**Files:**
- Create: `7.llm_dl_performance_course/modules/module01_models_and_hardware.md`
- Create: `7.llm_dl_performance_course/modules/module02_batching_and_memory.md`
- Create: `7.llm_dl_performance_course/modules/module03_concurrency_and_gpu_saturation.md`
- Create: `7.llm_dl_performance_course/modules/module04_sharding_and_mig.md`
- Create: `7.llm_dl_performance_course/modules/module05_parallelization.md`

**Step 1:** For each module file: title; 2–4 learning objectives; 1–2 paragraph summary; “Key concepts” list; “Curated resources” with explicit links (CS336, GPU MODE, NVIDIA, blogs from design); “Connection to biology (Part 2)” one short paragraph.

**Step 2:** Use exact resource URLs from user list where possible (buff.ly links can be kept or replaced with canonical URLs if known).

**Step 3:** Commit: `git add 7.llm_dl_performance_course/modules/ && git commit -m "feat(course7): add Part 1 core inference modules"`

---

## Task 4: Create Part 2 modules (biology application track)

**Files:**
- Create: `7.llm_dl_performance_course/modules/module06_single_cell_omics_dl.md`
- Create: `7.llm_dl_performance_course/modules/module07_inference_in_biology.md`

**Step 1:** `module06_single_cell_omics_dl.md`: gene/cell models (e.g. scGPT, Geneformer); where inference and perf matter in single-cell pipelines; link to Course 1–3 concepts (count matrices, clustering).

**Step 2:** `module07_inference_in_biology.md`: applying batching, memory, GPU saturation, MIG/sharding to biological inference; scaling and serving; large atlases.

**Step 3:** Commit: `git add 7.llm_dl_performance_course/modules/ && git commit -m "feat(course7): add Part 2 biology application modules"`

---

## Task 5: Curated readings and resources

**Files:**
- Create: `7.llm_dl_performance_course/resources/reading/curated_readings.md`

**Step 1:** Add two sections: “Courses / Videos” and “Repositories / Blogs”. Populate with user’s list:
- CS336 Stanford (Lecture 5–11): https://buff.ly/9o4jRxl
- GPU MODE: https://buff.ly/bPEHYL0
- LLM Inference Patterns: https://buff.ly/4ASDS5c
- LLM Inference (NVIDIA): https://buff.ly/muDMcor
- GPU Performance Engineering (E. Andere): https://buff.ly/sfGJNMK
- AI Performance Engineering (C. Fregly): https://buff.ly/8ey3GPa
- Inside NVIDIA GPUs (A. Gordic): https://buff.ly/wmgNwMg
- Understanding LLM Inference (article): https://buff.ly/NrNU3Qz

**Step 2:** Add a short note per resource and map to Part 1 / Part 2 (e.g. “Part 1 – batching, memory”).

**Step 3:** Commit: `git add 7.llm_dl_performance_course/resources/ && git commit -m "feat(course7): add curated readings"`

---

## Task 6: Quizzes and assignments placeholders

**Files:**
- Create: `7.llm_dl_performance_course/quizzes/quiz_module01.md` (Part 1 recap)
- Create: `7.llm_dl_performance_course/quizzes/quiz_part2.md` (Part 2 recap)
- Create: `7.llm_dl_performance_course/assignments/assignment1_inference_concepts.md`
- Create: `7.llm_dl_performance_course/assignments/assignment2_biology_application.md`

**Step 1:** Quizzes: 5–8 questions each (MC + short answer), 70% pass note; content aligned to module objectives.

**Step 2:** Assignments: brief prompt (e.g. “Explain batching and memory in your own words and give one biology example”; “Design a small deployment scenario for a cell-embedding model”).

**Step 3:** Commit: `git add 7.llm_dl_performance_course/quizzes/ 7.llm_dl_performance_course/assignments/ && git commit -m "feat(course7): add quizzes and assignment placeholders"`

---

## Task 7: Update root README and docs

**Files:**
- Modify: `README.md` — add Course 7 section
- Modify: `docs/MASTER_TRAINING_GUIDE.md` — add Course 7 to structure and prerequisites
- Modify: `config/notebooklm_courses.yaml` — add 7.llm_dl_performance_course entry

**Step 1:** In README.md, after Course 6 section, add “Course 7: LLMs, Deep Learning & Performance” with duration (e.g. 3–4 weeks), location `7.llm_dl_performance_course/`, START_HERE link, one-sentence description, “Part 1: Core inference / Part 2: Biology application”.

**Step 2:** In MASTER_TRAINING_GUIDE.md, add Course 7 to curriculum structure and prerequisite graph (after Course 1; 2/3 recommended).

**Step 3:** In config/notebooklm_courses.yaml, add key `7.llm_dl_performance_course` with name and sources (START_HERE, course_llm_dl_performance.md, modules/*.md, resources/reading/*.md, quizzes/*.md).

**Step 4:** Commit: `git add README.md docs/MASTER_TRAINING_GUIDE.md config/notebooklm_courses.yaml && git commit -m "feat(course7): integrate Course 7 into README, master guide, NotebookLM config"`

---

## Task 8: Optional labs and docs

**Files:**
- Create: `7.llm_dl_performance_course/labs/README.md` (placeholder: “Labs TBD; use curated resources and modules for now.”)
- Modify: `docs/ASSIGNMENT_RUBRICS.md` — add Course 7 row/section if file has per-course rubrics

**Step 1:** Add labs/README.md so START_HERE “labs” link is valid.

**Step 2:** If ASSIGNMENT_RUBRICS.md lists courses by number, add “Course 7” with placeholder criteria or “See assignment files.”

**Step 3:** Commit: `git add 7.llm_dl_performance_course/labs/README.md docs/ASSIGNMENT_RUBRICS.md && git commit -m "chore(course7): add labs placeholder and rubric placeholder"`

---

## Execution summary

- Tasks 1–8 in order; each task is 2–5 steps with exact paths.
- No automated tests; “verify” = file exists and links are consistent.
- After plan: offer Subagent-Driven (this session) vs Parallel Session (new session with executing-plans).
