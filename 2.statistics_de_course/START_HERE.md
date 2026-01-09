# START HERE: Statistics & Differential Expression

This is your **step-by-step path** through the course. Check items off as you complete them.

---

## 0) Orientation (5 minutes)

- [ ] Open `course_statistics_de.md` and skim the module list
- [ ] Open `resources/study_guide.md` to see the reference notes

---

## 1) Setup (30–90 minutes)

- [ ] Follow `setup.md` and confirm you can run R + Bioconductor packages
- [ ] Create a working folder (recommended):
  - [ ] `data/`
  - [ ] `results/`
  - [ ] `figures/`
- [ ] Copy templates:
  - [ ] `templates/design_matrix_template.md` → `results/design_notes.md`
  - [ ] `templates/de_report_template.md` → `results/de_report.md`

---

## 2) 4-Week Calendar (Follow This Exactly)

### Week 1 (Counts + Distributions + Design)

- [ ] **Day 1 (30–60m)**: Setup + create `data/`, `results/`, `figures/`
- [ ] **Day 2 (60–90m)**: **Lab 1** — `labs/lab01_explore_count_matrix.Rmd`
  - [ ] Save: dataset summary table + sample QC notes
- [ ] **Day 3 (60–90m)**: **Lab 2** — `labs/lab02_distributions.Rmd`
  - [ ] Save: mean–variance plots + Poisson vs NB notes
- [ ] **Day 4 (60–90m)**: **Lab 3** — `labs/lab03_design_matrices.Rmd`
  - [ ] Save: 3 example design matrices + contrast definitions to `results/design_notes.md`
- [ ] **Weekend (1–2h)**: **Assignment 1** — `assignments/assignment1_design_matrices.md`

### Week 2 (Normalization + Testing + Multiple Testing)

- [ ] **Day 1 (60–90m)**: **Lab 4** — `labs/lab04_normalization.Rmd`
  - [ ] Save: size factors summary + before/after normalization plots
- [ ] **Day 2 (60–90m)**: **Lab 5** — `labs/lab05_hypothesis_testing.Rmd`
- [ ] **Day 3 (45–75m)**: **Lab 6** — `labs/lab06_multiple_testing.Rmd`
  - [ ] Save: FDR simulation figure + threshold recommendation
- [ ] **Weekend (1–2h)**: **Assignment 2** — `assignments/assignment2_multiple_testing.md`

### Week 3 (DESeq2 + Effect Size + Tool Comparison)

- [ ] **Day 1 (45–75m)**: **Lab 7** — `labs/lab07_effect_size.Rmd`
  - [ ] Save: effect size vs significance examples
- [ ] **Day 2 (2–3h)**: **Lab 8** — `labs/lab08_deseq2.Rmd`
  - [ ] Save: results table (CSV) + MA plot + volcano plot
- [ ] **Day 3 (90–120m)**: **Lab 9** — `labs/lab09_compare_tools.Rmd`
  - [ ] Save: DESeq2 vs edgeR vs limma overlap table
- [ ] **Weekend (2–3h)**: **Assignment 3** — `assignments/assignment3_deseq2_analysis.md`

### Week 4 (scRNA-seq DE + Pseudobulk + Reporting)

- [ ] **Day 1 (60–90m)**: **Lab 10** — `labs/lab10_scrna_challenges.Rmd`
  - [ ] Save: pseudoreplication “bad example” notes
- [ ] **Day 2 (2–3h)**: **Lab 11** — `labs/lab11_pseudobulk.Rmd`
  - [ ] Save: pseudobulk DE results table (CSV)
- [ ] **Day 3 (60–90m)**: **Lab 12** — `labs/lab12_cell_level_de.Rmd`
  - [ ] Save: cell-level vs pseudobulk comparison plots
- [ ] **Day 4 (60–90m)**: **Lab 13** — `labs/lab13_visualization.Rmd`
  - [ ] Save: final figures for report
- [ ] **Weekend (2–3h)**: **Assignment 4** — `assignments/assignment4_pseudobulk_comparison.md`
- [ ] **Optional**: start **Final Project** — `assignments/final_project.md`

---

## 3) Labs (Reference List)

- [ ] Lab 1: `labs/lab01_explore_count_matrix.Rmd`
- [ ] Lab 2: `labs/lab02_distributions.Rmd`
- [ ] Lab 3: `labs/lab03_design_matrices.Rmd`
- [ ] Lab 4: `labs/lab04_normalization.Rmd`
- [ ] Lab 5: `labs/lab05_hypothesis_testing.Rmd`
- [ ] Lab 6: `labs/lab06_multiple_testing.Rmd`
- [ ] Lab 7: `labs/lab07_effect_size.Rmd`
- [ ] Lab 8: `labs/lab08_deseq2.Rmd`
- [ ] Lab 9: `labs/lab09_compare_tools.Rmd`
- [ ] Lab 10: `labs/lab10_scrna_challenges.Rmd`
- [ ] Lab 11: `labs/lab11_pseudobulk.Rmd`
- [ ] Lab 12: `labs/lab12_cell_level_de.Rmd`
- [ ] Lab 13: `labs/lab13_visualization.Rmd`

Suggested outputs to save in `results/`:
- [ ] session info (`sessionInfo()`)
- [ ] normalized counts summary
- [ ] MA plot + volcano plot + top DE table (CSV)

---

## 4) Quizzes (Quick Checks)

- [ ] `quizzes/quiz_module02.md`
- [ ] `quizzes/quiz_module05.md`
- [ ] `quizzes/quiz_module08.md`
- [ ] `quizzes/quiz_module11.md`

Target: **≥70%** before proceeding.

---

## 5) Assignments (Do After Relevant Labs)

- [ ] Assignment 1: `assignments/assignment1_design_matrices.md`
- [ ] Assignment 2: `assignments/assignment2_multiple_testing.md`
- [ ] Assignment 3: `assignments/assignment3_deseq2_analysis.md`
- [ ] Assignment 4: `assignments/assignment4_pseudobulk_comparison.md`

---

## 6) Final Project

- [ ] `assignments/final_project.md`
- [ ] Deliverables (minimum):
  - [ ] design matrix + contrast definitions
  - [ ] DE results table with FDR + log2FC
  - [ ] MA + volcano + heatmap
  - [ ] interpretation paragraph tying results to biology

---

## 7) Videos (Optional)

- [ ] Follow `resources/video_guide.md`

---

## If You Get Stuck

- [ ] Re-read “What Bad DE Analysis Looks Like” in `course_statistics_de.md`
- [ ] Confirm your replicate unit (sample/donor), and whether batch is confounded


