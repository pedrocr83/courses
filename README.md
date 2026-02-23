# Single-Cell RNA-seq Analysis Training Curriculum

**Comprehensive training program: From raw sequencing data to biological insight**

---

## 🎯 Program Overview

This curriculum provides complete training in single-cell RNA-seq analysis through **7 progressive courses** covering data processing, statistical analysis, clustering, integration, trajectory inference, cell-cell communication, and **LLMs/deep learning performance** (Course 7: intersection with biology).

**Total Duration:** 16-20 weeks (90-110 hours)  
**Format:** Self-paced with hands-on labs  
**Level:** Beginner → Advanced  
**Prerequisites:** Basic programming (R or Python), command line, genomics concepts

---

## 📚 Course Structure

### **Course 0: Prerequisites Assessment** (NEW - START HERE!)
**Duration:** Self-paced review  
**Location:** `0.prerequisites_assessment/`

Assess your readiness before starting:
- Programming skills (Python/R)
- Command line proficiency
- Biology & statistics fundamentals
- Git basics

**📍 START:** [`0.prerequisites_assessment/START_HERE.md`](0.prerequisites_assessment/START_HERE.md)

---

### **Course 1: Single-Cell RNA-seq Processing**
**Duration:** 3 weeks (15-20 hours)  
**Location:** `1.single_cell_processing_course/`

From raw FASTQ files to high-quality count matrices:
- Cell Ranger & kallisto|bustools pipelines
- Quality control and filtering
- Doublet detection
- Data export and reproducibility

**📍 START:** [`1.single_cell_processing_course/START_HERE.md`](1.single_cell_processing_course/START_HERE.md)

**🆕 NEW RESOURCES:**
- [QC Decision Framework](1.single_cell_processing_course/resources/qc_decision_framework.md) - Systematic QC threshold selection
- [Troubleshooting Guide](1.single_cell_processing_course/resources/troubleshooting_guide.md) - 60+ common errors solved
- [Self-Assessment Checklist](1.single_cell_processing_course/resources/self_assessment_checklist.md)

---

### **Course 2: Statistics & Differential Expression**
**Duration:** 4 weeks (20-25 hours)  
**Location:** `2.statistics_de_course/`

From count data to biological insight:
- Probability distributions for RNA-seq
- DESeq2 for bulk RNA-seq
- Pseudobulk analysis for scRNA-seq
- Multiple testing correction
- Effect size interpretation

**📍 START:** [`2.statistics_de_course/START_HERE.md`](2.statistics_de_course/START_HERE.md)

**🆕 NEW RESOURCES:**
- [Pseudobulk Theory Module](2.statistics_de_course/modules/module11_pseudobulk_theory.md) - Comprehensive treatment of pseudoreplication

---

### **Course 3: Clustering & Cell Type Annotation**
**Duration:** 3 weeks (15-20 hours)  
**Location:** `3.clustering_annotation_course/`

From processed counts to biological identity:
- Dimensionality reduction (PCA, UMAP, t-SNE)
- Graph-based clustering (Leiden/Louvain)
- Marker gene identification
- Manual and automated annotation
- Annotation confidence assessment

**📍 START:** [`3.clustering_annotation_course/START_HERE.md`](3.clustering_annotation_course/START_HERE.md)

**🆕 NEW RESOURCES:**
- [Annotation Confidence Framework](3.clustering_annotation_course/resources/annotation_confidence_framework.md) - 3-tier confidence system with evidence scoring
- [Self-Assessment Checklist](3.clustering_annotation_course/resources/self_assessment_checklist.md)

---

### **Course 4: Data Integration & Batch Correction**
**Duration:** 3 weeks (15-20 hours)  
**Location:** `4.data_integration_course/`

Combining multiple datasets into unified atlases:
- Batch effect diagnosis
- Integration methods (Harmony, Seurat, scVI)
- Integration quality evaluation
- Reference mapping

**📍 START:** [`4.data_integration_course/START_HERE.md`](4.data_integration_course/START_HERE.md)

**🆕 NEW RESOURCES:**
- [Integration Method Decision Guide](4.data_integration_course/resources/integration_method_decision_guide.md) - Interactive decision tree + detailed comparisons

---

### **Course 5: Trajectory & Pseudotime Analysis**
**Duration:** 3 weeks (15-20 hours)  
**Location:** `5.trajectory_analysis_course/`

Understanding cellular dynamics and differentiation:
- Pseudotime inference (DPT, PAGA, Monocle, Slingshot)
- RNA velocity (scVelo)
- Fate prediction (CellRank)
- Trajectory differential expression
- **Validation strategies**

**📍 START:** [`5.trajectory_analysis_course/START_HERE.md`](5.trajectory_analysis_course/START_HERE.md)

**🆕 NEW RESOURCES:**
- [Trajectory Validation Module](5.trajectory_analysis_course/modules/module11_trajectory_validation.md) - Multi-level validation framework

---

### **Course 6: Cell-Cell Communication Analysis**
**Duration:** 4 weeks (20-25 hours)  
**Location:** `6.cell_cell_communication_course/`

From expression to interaction networks:
- CCC inference tools (CellPhoneDB, CellChat, NicheNet)
- Method comparison and consensus (LIANA)
- Network analysis
- Spatial context integration
- **Experimental validation design**

**📍 START:** [`6.cell_cell_communication_course/START_HERE.md`](6.cell_cell_communication_course/START_HERE.md)

**🆕 NEW RESOURCES:**
- [CCC Validation Module](6.cell_cell_communication_course/modules/module13_ccc_validation.md) - From prediction to experimental validation

---

### **Course 7: LLMs, Deep Learning & Performance**
**Duration:** 3–4 weeks (15–22 hours)  
**Location:** `7.llm_dl_performance_course/`

Efficient inference and AI performance at the intersection with biology:
- **Part 1:** Core inference (batching, memory, concurrency, GPU saturation, sharding/MIG, parallelization)
- **Part 2:** Biology application track (single-cell & omics DL; applying inference to gene/cell models)

**Prerequisites:** Course 1; Courses 2/3 recommended.

**📍 START:** [`7.llm_dl_performance_course/START_HERE.md`](7.llm_dl_performance_course/START_HERE.md)

**Resources:**
- [Curated readings](7.llm_dl_performance_course/resources/reading/curated_readings.md) — courses, videos, blogs (CS336, GPU MODE, NVIDIA, inference patterns, perf engineering)

---

## 🆕 What's New: Phase 1 Improvements

### Major Enhancements

1. **Prerequisites Assessment (Course 0)** - Assess readiness before starting
2. **Decision Frameworks** - Systematic guidance for QC, integration, and annotation choices
3. **Validation Modules** - Critical evaluation of trajectories and CCC predictions
4. **Expanded Pseudobulk** - Comprehensive treatment of pseudoreplication
5. **Assignment Rubrics** - Transparent grading criteria for all 20 assignments
6. **Troubleshooting Guide** - Self-service problem solving (60+ errors)
7. **Self-Assessment Checklists** - Track your mastery systematically

### Key Documents

- **[Master Training Guide](docs/MASTER_TRAINING_GUIDE.md)** - Complete curriculum roadmap
- **[Assignment Rubrics](docs/ASSIGNMENT_RUBRICS.md)** - Grading criteria for all assignments
- **[Improvement Suggestions](docs/COURSE_IMPROVEMENT_SUGGESTIONS.md)** - Roadmap for future enhancements
- **[Implementation Summary](docs/IMPLEMENTATION_SUMMARY.md)** - What we've built and why

---

## 🚀 Quick Start

### For New Students

**Step 1: Assess Prerequisites (Week 0)**
```
1. Review 0.prerequisites_assessment/START_HERE.md
2. Take self-assessments
3. Build missing skills if needed
4. Set up computational environment
```

**Step 2: Follow Course Sequence**
```
Course 1 → Course 2 + Course 3 (can overlap) → Course 4 → Course 5 → Course 6
```

**Step 3: Use Support Materials**
```
- Decision frameworks BEFORE key choices
- Troubleshooting guides WHEN stuck
- Self-assessment checklists REGULARLY
- Assignment rubrics BEFORE submitting
```

### For Instructors

**Before Teaching:**
1. Review all decision frameworks
2. Familiarize with assignment rubrics
3. Note new modules (Pseudobulk, Validation × 2)
4. Identify which resources to emphasize

**During Teaching:**
- Point students to frameworks at decision points
- Use rubrics for consistent grading
- Monitor troubleshooting guide usage
- Collect feedback for Phase 2

---

## 📊 Learning Paths

### 4-Month Intensive (20-25 hrs/week)
- **Month 1:** Course 1 + Course 2
- **Month 2:** Course 3 + Course 4
- **Month 3:** Course 5
- **Month 4:** Course 6

### 6-Month Standard (12-15 hrs/week)
- **Months 1-2:** Course 1
- **Months 2-3:** Course 2
- **Months 3-4:** Course 3
- **Months 4-5:** Course 4
- **Months 5-6:** Course 5
- **Month 6:** Course 6

### 9-Month Relaxed (8-10 hrs/week)
- Linear progression with extra practice time

---

## 🎯 Learning Objectives

By completing this curriculum, you will be able to:

**Technical Skills:**
- ✅ Process raw scRNA-seq data end-to-end
- ✅ Perform rigorous differential expression analysis
- ✅ Cluster cells and annotate cell types
- ✅ Integrate multi-batch datasets
- ✅ Infer developmental trajectories
- ✅ Analyze cell-cell communication networks

**Conceptual Understanding:**
- ✅ Make justified analysis decisions
- ✅ Identify common pitfalls and artifacts
- ✅ Critically evaluate methods
- ✅ Validate computational predictions
- ✅ Report results with appropriate confidence

**Professional Skills:**
- ✅ Build reproducible analysis pipelines
- ✅ Create publication-quality figures
- ✅ Document analyses thoroughly
- ✅ Design validation experiments
- ✅ Communicate findings effectively

---

## 📁 Repository Structure

```
courses/
├── 0.prerequisites_assessment/          [NEW]
│   ├── START_HERE.md
│   └── assessments/
├── 1.single_cell_processing_course/
│   ├── START_HERE.md
│   ├── course_single_cell_processing.md
│   ├── labs/
│   ├── assignments/
│   ├── quizzes/
│   └── resources/                       [ENHANCED]
│       ├── qc_decision_framework.md     [NEW]
│       ├── troubleshooting_guide.md     [NEW]
│       └── self_assessment_checklist.md [NEW]
├── 2.statistics_de_course/
│   └── modules/
│       └── module11_pseudobulk_theory.md [NEW]
├── 3.clustering_annotation_course/
│   └── resources/
│       ├── annotation_confidence_framework.md [NEW]
│       └── self_assessment_checklist.md [NEW]
├── 4.data_integration_course/
│   └── resources/
│       └── integration_method_decision_guide.md [NEW]
├── 5.trajectory_analysis_course/
│   └── modules/
│       └── module11_trajectory_validation.md [NEW]
├── 6.cell_cell_communication_course/
│   └── modules/
│       └── module13_ccc_validation.md   [NEW]
├── 7.llm_dl_performance_course/          [NEW]
│   ├── START_HERE.md
│   ├── course_llm_dl_performance.md
│   ├── setup.md
│   ├── glossary.md
│   ├── modules/                          (Part 1: core inference; Part 2: biology)
│   ├── resources/reading/curated_readings.md
│   ├── quizzes/
│   ├── assignments/
│   └── labs/README.md
├── config/                               [NEW]
│   ├── notebooklm_courses.yaml          # NotebookLM upload config
│   └── quiz_template_detailed.md        # Template for detailed quizzes
├── scripts/                              [NEW]
│   ├── notebooklm_upload.py            # Upload courses to NotebookLM
│   └── requirements-notebooklm.txt     # Optional: NotebookLM deps
├── [course]/resources/reading/           [NEW - per course]
│   └── curated_readings.md             # Recommended papers
└── docs/
    ├── MASTER_TRAINING_GUIDE.md
    ├── COURSE_IMPROVEMENT_SUGGESTIONS.md
    ├── ASSIGNMENT_RUBRICS.md            [NEW]
    ├── IMPROVEMENTS_IMPLEMENTED.md      [NEW]
    ├── STANDARDIZED_START_HERE_TEMPLATE.md [NEW]
    └── IMPLEMENTATION_SUMMARY.md        [NEW]
```

---

## 🛠️ Software Requirements

### R Environment
```r
# Core packages
BiocManager::install(c("DESeq2", "edgeR", "limma", "SingleR", 
                       "Seurat", "batchelor", "SingleCellExperiment"))
```

### Python Environment
```bash
pip install scanpy anndata scvi-tools harmonypy cellrank scvelo \
            leidenalg liana
```

### Pipeline Tools
- Cell Ranger 8.x
- STAR 2.7.x
- kallisto + bustools

**Detailed setup:** See `setup.md` in each course folder

---

## 🎓 Assessments & Certification

### Assessment Structure
- **Quizzes:** Module-level knowledge checks (≥70% to pass)
- **Assignments:** 20 practical assignments across 6 courses
- **Final Projects:** End-of-course capstone projects (6 total)
- **Capstone:** Integrated analysis using all skills

### Certification Requirements
- [ ] Complete all 6 courses
- [ ] Pass all quizzes (≥70%)
- [ ] Submit all assignments
- [ ] Complete all final projects
- [ ] Capstone project (optional but recommended)

---

## 📖 Key Resources

### Decision Frameworks (Use Before Making Choices)
- [QC Decision Framework](1.single_cell_processing_course/resources/qc_decision_framework.md)
- [Integration Method Decision Guide](4.data_integration_course/resources/integration_method_decision_guide.md)
- [Annotation Confidence Framework](3.clustering_annotation_course/resources/annotation_confidence_framework.md)

### Validation Modules (Use for Advanced Analyses)
- [Trajectory Validation](5.trajectory_analysis_course/modules/module11_trajectory_validation.md)
- [CCC Validation](6.cell_cell_communication_course/modules/module13_ccc_validation.md)

### Support Materials
- [Troubleshooting Guide (Course 1)](1.single_cell_processing_course/resources/troubleshooting_guide.md)
- [Assignment Rubrics (All Courses)](docs/ASSIGNMENT_RUBRICS.md)
- [Self-Assessment Checklists](1.single_cell_processing_course/resources/self_assessment_checklist.md)

### Comprehensive Guides
- [Master Training Guide](docs/MASTER_TRAINING_GUIDE.md) - Complete program overview
- [Improvement Suggestions](docs/COURSE_IMPROVEMENT_SUGGESTIONS.md) - Future enhancements

---

## 🎧 Listen & Learn with NotebookLM

Upload course content to [Google NotebookLM](https://notebooklm.google.com) to generate **audio overviews (podcasts)**, quizzes, and flashcards. Listen while commuting or ask questions about your materials.

### Setup

```bash
pip install -r scripts/requirements-notebooklm.txt
playwright install chromium
notebooklm login   # one-time Google auth (opens browser)
```

### Usage

```bash
# Upload a single course
python3 scripts/notebooklm_upload.py 1.single_cell_processing_course

# Upload all courses
python3 scripts/notebooklm_upload.py all

# Upload and generate a podcast
python3 scripts/notebooklm_upload.py 1.single_cell_processing_course --generate audio

# Preview files without uploading
python3 scripts/notebooklm_upload.py 1.single_cell_processing_course --dry-run
```

### What Gets Uploaded

Each course uploads: START_HERE, study guides, modules, resources, **reading material** (`resources/reading/`), and **quizzes**. Edit `config/notebooklm_courses.yaml` to customize sources.

### Reading Material

Each course has `resources/reading/curated_readings.md` with recommended papers. Add more `.md` files or PDFs to `[course]/resources/reading/`—they will be included when you upload.

### Detailed Quizzes

All module quizzes use a detailed format (15 questions: 10 MC + 5 short answer). Template: `config/quiz_template_detailed.md`.

### Alternative: Claude Skill

[anything-to-notebooklm](https://github.com/joeseesun/anything-to-notebooklm) – Claude Code skill for turning URLs, PDFs, and markdown into podcasts/quizzes via natural language.

---

## 🤝 Getting Help

### Self-Service Resources
1. **Troubleshooting guide** - Check for your error
2. **Decision frameworks** - Guidance for key choices
3. **Self-assessment checklists** - Identify knowledge gaps
4. **Course discussion forums** - Search previous questions

### When to Ask for Help
- After checking troubleshooting guide
- After reviewing relevant course materials
- When you have a clear, specific question
- With error messages and what you've tried

### How to Ask Good Questions
Include:
- Clear problem description
- Full error message
- Software versions
- What you've already tried
- Minimal code to reproduce

---

## 🔄 Continuous Improvement

This curriculum is actively maintained and improved. **Phase 1 (HIGH PRIORITY)** improvements complete. **Phase 2** in planning:

### Coming Next (Phase 2)
- Video walkthroughs for all modules
- Integrated case study dataset (one dataset, six courses)
- Additional troubleshooting guides (Courses 2-6)
- More self-assessment checklists
- Cloud computing options

### Contribute
- Report errors or unclear sections
- Suggest improvements
- Share successful teaching strategies
- Contribute code examples

**See:** [`docs/COURSE_IMPROVEMENT_SUGGESTIONS.md`](docs/COURSE_IMPROVEMENT_SUGGESTIONS.md)

---

## 📊 Success Metrics

**Student Outcomes:**
- Course completion rate: Target 80%
- Quiz performance: Target ≥75% average
- Assignment quality: Target 85% proficient level
- Time to competency: 4-9 months depending on pace

**Learning Impact:**
- Can independently analyze scRNA-seq data
- Produces publication-quality work
- Makes justified, documented decisions
- Appropriately validates findings

---

## 📜 Citation

If you use these materials in teaching or research, please cite:

```
[Citation to be added - course repository or publication]
```

---

## 📞 Contact & Support

**Course Administrators:** [Contact info]  
**Technical Issues:** [GitHub issues link]  
**General Questions:** [Discussion forum]  
**Contributions:** [Pull requests welcome]

---

## 📄 License

[Specify license - e.g., CC BY-NC-SA 4.0 for educational materials]

---

## 🙏 Acknowledgments

**Development Team:**
- [List contributors]

**Reviewers & Testers:**
- [List pilot students and reviewers]

**Funding:**
- [Funding sources if applicable]

**Tools & Resources:**
- Scanpy, Seurat, DESeq2, and all open-source tools used
- Data providers (10x Genomics, public repositories)
- Broader single-cell community

---

## 🌟 What Makes This Curriculum Unique

1. **Comprehensive Coverage** - All analysis stages, not just parts
2. **Decision Frameworks** - Systematic guidance, not arbitrary choices
3. **Validation Focus** - Critical evaluation, not just prediction
4. **Statistical Rigor** - Proper methods (e.g., pseudobulk), not just convenience
5. **Transparent Assessment** - Clear rubrics, not subjective grading
6. **Self-Guided Learning** - Resources for independent progress
7. **Publication Quality** - Standards for real research, not just exercises

---

**Ready to begin your single-cell journey?**

**👉 START HERE:** [`0.prerequisites_assessment/START_HERE.md`](0.prerequisites_assessment/START_HERE.md)

---

**Last Updated:** February 2, 2026  
**Version:** 2.0 (Phase 1 Improvements Complete)  
**Status:** ✅ Ready for Use
