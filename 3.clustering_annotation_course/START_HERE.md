# START HERE: Clustering & Cell Type Annotation

This is your **step-by-step path** through the course. Check items off as you complete them.

---

## 1) Setup

- [ ] Follow `setup.md` (Python Scanpy recommended; R optional)
- [ ] Create folders:
  - [ ] `data/`
  - [ ] `results/`
  - [ ] `figures/`
- [ ] Copy the report template:
  - [ ] `templates/annotation_report_template.md` → `results/annotation_report.md`

---

## 2) 3-Week Calendar (Follow This Exactly)

### Week 1 (Dimensionality Reduction)

- [ ] **Day 1 (30–60m)**: Setup + create `data/`, `results/`, `figures/`
- [ ] **Day 2 (45–60m)**: **Lab 1** — `labs/lab01_high_dim_explore.ipynb`
  - [ ] Save: `results/lab01_highdim_notes.md`
- [ ] **Day 3 (60–90m)**: **Lab 2** — `labs/lab02_pca.ipynb`
  - [ ] Save: elbow plot + chosen PCs
- [ ] **Day 4 (60–90m)**: **Lab 3** — `labs/lab03_umap_tsne.ipynb`
  - [ ] Save: UMAP/t-SNE parameter comparison figures
- [ ] **Weekend (1–2h)**: Assignment 1

### Week 2 (Clustering)

- [ ] **Day 1 (45–60m)**: **Lab 4** — `labs/lab04_cell_cell_distances.ipynb`
  - [ ] Save: `results/lab04_distance_notes.md`
- [ ] **Day 2 (45–75m)**: **Lab 5A** — `labs/lab05A_build_graph_cluster.ipynb`
  - [ ] Save: `results/lab05A_clusters.csv`
- [ ] **Day 3 (45–75m)**: **Lab 5B** — `labs/lab05B_resolution_effects.ipynb`
  - [ ] Save: `figures/lab05B_resolution_grid.png`
- [ ] **Day 4 (45–60m)**: **Lab 6** — `labs/lab06_cluster_evaluation.ipynb`
  - [ ] Save: `results/lab06_cluster_eval.md`
- [ ] **Weekend (1–2h)**: Assignment 2

### Week 3 (Annotation)

- [ ] **Day 1 (60–90m)**: **Lab 7** — `labs/lab07_markers.ipynb`
  - [ ] Save: marker dotplots/heatmaps
- [ ] **Day 2 (60–90m)**: **Lab 8** — `labs/lab08_manual_annotation.ipynb`
  - [ ] Save: `results/lab08_manual_annotations.csv` + `figures/lab08_umap_celltypes.png`
- [ ] **Day 3 (60–90m)**: **Lab 9** — `labs/lab09_annotation.ipynb`
  - [ ] Save: automated labels + confidence summaries
- [ ] **Day 4 (45–75m)**: **Lab 10** — `labs/lab10_refinement_qc.ipynb`
  - [ ] Save: `results/lab10_final_labels.csv`
- [ ] **Weekend (2–3h)**: Assignment 3 + start Final Project

---

## 3) Labs (Reference List)

- [ ] `labs/lab01_high_dim_explore.ipynb`
- [ ] `labs/lab02_pca.ipynb`
- [ ] `labs/lab03_umap_tsne.ipynb`
- [ ] `labs/lab04_cell_cell_distances.ipynb`
- [ ] `labs/lab05A_build_graph_cluster.ipynb`
- [ ] `labs/lab05B_resolution_effects.ipynb`
- [ ] `labs/lab06_cluster_evaluation.ipynb`
- [ ] `labs/lab07_markers.ipynb`
- [ ] `labs/lab08_manual_annotation.ipynb`
- [ ] `labs/lab09_annotation.ipynb`
- [ ] `labs/lab10_refinement_qc.ipynb`

---

## 4) Quizzes

- [ ] `quizzes/quiz_module02.md` (PCA)
- [ ] `quizzes/quiz_module05.md` (Graph clustering)
- [ ] `quizzes/quiz_module08.md` (Manual annotation)

---

## 5) Assignments

- [ ] `assignments/assignment1_dimensionality_reduction.md`
- [ ] `assignments/assignment2_clustering.md`
- [ ] `assignments/assignment3_annotation.md`

---

## 6) Final Project

- [ ] `assignments/final_project.md`
- [ ] Deliverables:
  - [ ] `results/annotation_report.md` completed (with evidence + confidence)
  - [ ] figures folder with final plots
  - [ ] exported `annotations.csv`

---

## If You Get Stuck

- [ ] Read “What Bad Clustering / Bad Annotation Looks Like” in `course_clustering_annotation.md`
- [ ] Check: are clusters driven by QC metrics? do clusters have specific markers?


