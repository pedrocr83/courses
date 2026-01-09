# START HERE: Data Integration & Batch Correction

This is your **step-by-step path** through the course. Check items off as you complete them.

---

## 1) Setup

- [ ] Follow `setup.md`
- [ ] Create folders:
  - [ ] `data/`
  - [ ] `results/`
  - [ ] `figures/`
- [ ] Copy the report template:
  - [ ] `templates/integration_report_template.md` → `results/integration_report.md`

---

## 2) Labs (Do in Order)

- [ ] `labs/lab02_batch_diagnosis.ipynb` (visual + LISI-style metrics)
- [ ] `labs/lab05_harmony.ipynb`
- [ ] `labs/lab07_scvi.ipynb` (optional if scvi-tools installed)
- [ ] `labs/lab08_evaluation.ipynb`
- [ ] `labs/lab10_reference_mapping.ipynb`

Save these outputs:
- [ ] pre-integration UMAP by batch + by cell type
- [ ] post-integration UMAP by batch + by cell type
- [ ] table of metrics (mixing + conservation)
- [ ] runtime/parameter notes

---

## 3) Quizzes

- [ ] `quizzes/quiz_module02.md`
- [ ] `quizzes/quiz_module05.md`
- [ ] `quizzes/quiz_module08.md`

---

## 4) Assignments

- [ ] `assignments/assignment1_batch_diagnosis.md`
- [ ] `assignments/assignment2_method_comparison.md`
- [ ] `assignments/assignment3_evaluation_report.md`

---

## 5) Final Project

- [ ] `assignments/final_project.md`
- [ ] Deliverables:
  - [ ] integrated object (`integrated_adata.h5ad` or Seurat object)
  - [ ] `results/integration_report.md` completed
  - [ ] method comparison table + figures

---

## If You Get Stuck

- [ ] Read “What Bad Integration Looks Like” in `course_data_integration.md`
- [ ] Check: are you over-correcting (biology lost) or under-correcting (batches remain)?


