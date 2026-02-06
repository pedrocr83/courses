# START HERE: Trajectory & Pseudotime Analysis

This is your **step-by-step path** through the course. Check items off as you complete them.

> **Updated February 2026** — Now includes CellRank 2 multiview fate mapping (Week 3).

---

## 1) Setup

- [ ] Follow `setup.md` (Python + CellRank + scVelo recommended)
- [ ] Clone the CellRank Protocol repo for reference: `git clone https://github.com/theislab/cellrank_protocol.git`
- [ ] Create folders:
  - [ ] `data/`
  - [ ] `results/`
  - [ ] `figures/`
- [ ] Copy the report template:
  - [ ] `templates/trajectory_report_template.md` → `results/trajectory_report.md`
- [ ] Verify installation:

```python
import cellrank as cr
from cellrank.kernels import VelocityKernel, CytoTRACEKernel, PseudotimeKernel, RealTimeKernel
print(f"CellRank {cr.__version__} — all kernels available!")
```

---

## 2) Week 1: Foundations (Labs 1-3)

- [ ] Read Modules 1-3 in `course_trajectory_analysis.md`
- [ ] `labs/lab03_diffusion_pt.ipynb`
- [ ] `quizzes/quiz_module02.md`

**Save these outputs:**
- [ ] pseudotime UMAP + pseudotime-by-cluster summary

---

## 3) Week 2: Trajectory Inference Methods (Labs 4-6)

- [ ] Read Modules 4-6 in `course_trajectory_analysis.md`
- [ ] `labs/lab06_velocity.ipynb`
- [ ] `quizzes/quiz_module06.md`
- [ ] `assignments/assignment1_pseudotime.md`
- [ ] `assignments/assignment2_velocity.md`

**Save these outputs:**
- [ ] velocity stream plot + velocity gene plots
- [ ] comparison of DPT, Monocle3, and Slingshot pseudotime orderings

---

## 4) Week 3: CellRank 2 — Multiview Fate Mapping (Labs 8-9D)

> This is the new core of the course. Based on the [CellRank Nature Protocols paper](https://doi.org/10.1038/s41596-025-01314-w) (Weiler & Theis, 2026).

- [ ] Read Modules 7-12 in `course_trajectory_analysis.md`
- [ ] Read `resources/cellrank2_multiview_guide.md` for kernel reference
- [ ] **Lab 8**: `labs/lab09_cellrank.ipynb` — VelocityKernel fate mapping
- [ ] **Lab 9A**: `labs/lab09A_cellrank_cytotrace_kernel.ipynb` — CytoTRACEKernel
- [ ] **Lab 9B**: `labs/lab09B_cellrank_pseudotime_kernel.ipynb` — PseudotimeKernel
- [ ] **Lab 9C**: `labs/lab09C_cellrank_realtime_kernel.ipynb` — RealTimeKernel
- [ ] **Lab 9D**: `labs/lab09D_cellrank_kernel_combination.ipynb` — Kernel combination
- [ ] `quizzes/quiz_module07_08_cellrank2.md`
- [ ] `quizzes/quiz_module09_12_cellrank2.md`

**Save these outputs:**
- [ ] fate probability UMAPs from each kernel
- [ ] kernel comparison correlation matrix
- [ ] terminal states from combined kernel analysis
- [ ] driver gene lists per lineage

---

## 5) Week 4: Downstream Analysis & Comparison (Labs 10-12)

- [ ] Read Modules 13-15 in `course_trajectory_analysis.md`
- [ ] `labs/lab07_trajectory_de.ipynb`
- [ ] `labs/lab10_comparison.ipynb`
- [ ] `assignments/assignment3_trajectory_de.md`
- [ ] `assignments/assignment4_cellrank_multiview.md`

**Save these outputs:**
- [ ] top pseudotime-associated genes + trend plots
- [ ] CellRank driver genes + gene trend plots
- [ ] method comparison notes (what changes, what stays stable)

---

## 6) Final Project

- [ ] `assignments/final_project.md`
- [ ] Deliverables:
  - [ ] `results/trajectory_report.md` completed
  - [ ] pseudotime + velocity figures
  - [ ] CellRank 2 multiview analysis (≥2 kernels + combined)
  - [ ] fate probability figures + terminal state validation
  - [ ] driver gene lists (trajectory DE + CellRank driver genes)
  - [ ] kernel comparison & agreement analysis

---

## Key Resources

| Resource | Location |
|----------|----------|
| Course outline | `course_trajectory_analysis.md` |
| CellRank 2 guide | `resources/cellrank2_multiview_guide.md` |
| CellRank Protocol | [github.com/theislab/cellrank_protocol](https://github.com/theislab/cellrank_protocol) |
| Protocol Jupyter Book | [theislab.github.io/cellrank_protocol](https://theislab.github.io/cellrank_protocol/index.html) |
| Nature Protocols paper | [doi:10.1038/s41596-025-01314-w](https://doi.org/10.1038/s41596-025-01314-w) |
| CellRank docs | [cellrank.readthedocs.io](https://cellrank.readthedocs.io/) |
| Protocol data | [Figshare](https://doi.org/10.6084/m9.figshare.c.7752290.v1) |

---

## If You Get Stuck

- [ ] Read "What Bad Trajectory Analysis Looks Like" in `course_trajectory_analysis.md`
- [ ] Read the CellRank 2 kernel decision guide in `course_trajectory_analysis.md`
- [ ] Check the [CellRank discourse forum](https://discourse.scverse.org/c/ecosystem/cellrank/)
- [ ] Ask: is your system truly continuous, or are you forcing a trajectory on discrete types?
- [ ] Ask: have you tried a different kernel to see if the result is robust?
