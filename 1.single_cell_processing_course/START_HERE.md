# START HERE: Single-Cell RNA-seq Sample Processing

This is your **step-by-step path** through the course. Check items off as you complete them.

---

## 0) Folder Orientation (2 minutes)

- [ ] Open `course_single_cell_processing.md` and skim the **Course Overview** + **Course Structure**
- [ ] Open `resources/study_guide.md` to see the original reference notes

---

## 1) Setup (30–90 minutes)

- [ ] Follow `setup.md` and install required software (at minimum: Python + `scanpy`)
- [ ] Create a working folder (recommended):
  - [ ] `data/raw/`
  - [ ] `data/processed/`
  - [ ] `results/`
  - [ ] `figures/`
- [ ] Copy templates you will use:
  - [ ] `templates/metadata_template.csv` → your `data/metadata.csv`
  - [ ] `templates/processing_log_template.md` → your `results/processing_log.md`
  - [ ] `templates/qc_report_template.md` → your `results/qc_report.md`

---

## 2) 3-Week Calendar (Follow This Exactly)

### Week 1 (Foundations + FASTQ)

- [ ] **Day 1 (30–60m)**: Setup + folder creation + copy templates
- [ ] **Day 2 (45–60m)**: **Lab 1** — `labs/lab01_explore_dataset.ipynb`
  - [ ] Save: dataset summary (cells/genes), basic QC metrics table to `results/`
- [ ] **Day 3 (30–45m)**: **Lab 2** — `labs/lab02_metadata_template.ipynb`
  - [ ] Save: `data/metadata.csv` + `results/metadata_validation.md`
- [ ] **Day 4 (60–90m)**: **Lab 3** — `labs/lab03_fastq_inspection.ipynb`
  - [ ] Save: FASTQ structure notes + a short “read anatomy” section for your chemistry
- [ ] **Weekend (1–2h)**: **Assignment 1** — `assignments/assignment1_fastq_report.md`

### Week 2 (Reference + Quantification + Cell Calling)

- [ ] **Day 1 (45–60m)**: **Lab 4** — `labs/lab04_reference_prep.ipynb`
  - [ ] Save: `results/reference_notes.md` (URLs, versions, checksums)
- [ ] **Day 2 (2–6h, optional heavy compute)**: **Lab 5A** — `labs/lab05A_cellranger.ipynb`
  - [ ] Save: `results/lab05A_cellranger_notes.md` + copy metrics into `results/processing_log.md`
- [ ] **Day 3 (2–6h, optional heavy compute)**: **Lab 5B** — `labs/lab05B_kallisto_bustools.ipynb`
  - [ ] Save: `results/lab05B_kb_notes.md` + copy metrics into `results/processing_log.md`
- [ ] **Day 4 (60–90m)**: **Lab 6** — `labs/lab06_cell_calling.ipynb`
  - [ ] Save: knee plot + your chosen “cell calling” rule and rationale
- [ ] **Weekend (2–3h)**: **Assignment 2** — `assignments/assignment2_pipeline_comparison.md`

### Week 3 (Matrix + QC + Filtering + Export)

- [ ] **Day 1 (45–60m)**: **Lab 7** — `labs/lab07_count_matrix_loading.ipynb`
  - [ ] Save: matrix shape, sparsity, raw vs filtered notes
- [ ] **Day 2 (60–90m)**: **Lab 8** — `labs/lab08_qc_metrics.ipynb`
  - [ ] Save: QC distributions, proposed thresholds/outlier method
- [ ] **Day 3 (60–120m)**: **Lab 9** — `labs/lab09_filtering_doublets.ipynb`
  - [ ] Save: before/after metrics + doublet decision log
- [ ] **Day 4 (30–60m)**: **Lab 10** — `labs/lab10_export_reproducibility.ipynb`
  - [ ] Save: `results/repro_bundle/` with exports + versions + params
- [ ] **Weekend (2–3h)**: **Assignment 3** — `assignments/assignment3_qc_report.md`
- [ ] **Optional**: start **Final Project** — `assignments/final_project.md`

---

## 3) Labs (Reference List)

### Week 1
- [ ] **Lab 1**: `labs/lab01_explore_dataset.ipynb`
- [ ] **Lab 2**: `labs/lab02_metadata_template.ipynb`
- [ ] **Lab 3**: `labs/lab03_fastq_inspection.ipynb`

### Week 2
- [ ] **Lab 4**: `labs/lab04_reference_prep.ipynb`
- [ ] **Lab 5A**: `labs/lab05A_cellranger.ipynb`
- [ ] **Lab 5B**: `labs/lab05B_kallisto_bustools.ipynb`
- [ ] **Lab 6**: `labs/lab06_cell_calling.ipynb`

### Week 3
- [ ] **Lab 7**: `labs/lab07_count_matrix_loading.ipynb`
- [ ] **Lab 8**: `labs/lab08_qc_metrics.ipynb`
- [ ] **Lab 9**: `labs/lab09_filtering_doublets.ipynb`
- [ ] **Lab 10**: `labs/lab10_export_reproducibility.ipynb`

---

## 4) Quizzes (Quick Checks)

- [ ] `quizzes/quiz_module01.md`
- [ ] `quizzes/quiz_module05.md`
- [ ] `quizzes/quiz_module08.md`

Target: **≥70%** before proceeding.

---

## 5) Assignments (Do After Relevant Labs)

- [ ] **Assignment 1**: `assignments/assignment1_fastq_report.md`
  - Output: FASTQ inspection report (attach figures/snippets)
- [ ] **Assignment 2**: `assignments/assignment2_pipeline_comparison.md`
  - Output: comparison table + pros/cons + reproducibility notes
- [ ] **Assignment 3**: `assignments/assignment3_qc_report.md`
  - Output: completed `results/qc_report.md` + figures + threshold rationale

---

## 6) Final Project (Capstone)

- [ ] `assignments/final_project.md`
- [ ] Deliverables:
  - [ ] a filtered count matrix (e.g., `.h5ad`/`.mtx`)
  - [ ] QC report with **before/after** plots
  - [ ] processing log with versions + parameters

---

## 7) Watch Videos (Optional but Helpful)

- [ ] Follow `resources/video_guide.md` by module/week

---

## If You Get Stuck

- [ ] Re-read: **“What Bad Processing / Bad QC Looks Like”** section in `course_single_cell_processing.md`
- [ ] Keep notes in `results/processing_log.md` (what you tried, what changed, what you observed)


