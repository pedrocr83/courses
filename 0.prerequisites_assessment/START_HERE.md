# Prerequisites Assessment

**Complete this module before starting Course 1 to ensure you're ready for the training program.**

---

## What This Module Covers

This self-paced assessment helps you:
1. Evaluate your readiness for the single-cell RNA-seq curriculum
2. Identify skill gaps that need attention
3. Access resources to build missing prerequisites
4. Set realistic expectations for time commitment

**Time Required:** 3-5 hours (assessment + skill building as needed)

---

## Assessment Areas

### 1. Programming Skills
- [ ] Complete `assessments/programming_assessment.md`
- [ ] Minimum passing: 70% (7/10 questions correct)
- [ ] If you score < 70%: Review `resources/programming_primer.md`

### 2. Command Line Basics
- [ ] Complete `assessments/command_line_assessment.md`
- [ ] Minimum passing: 70% (7/10 questions correct)
- [ ] If you score < 70%: Review `resources/command_line_primer.md`

### 3. Biology & Statistics Fundamentals
- [ ] Complete `assessments/biology_stats_assessment.md`
- [ ] Minimum passing: 70% (14/20 questions correct)
- [ ] If you score < 70%: Review `resources/biology_stats_primer.md`

### 4. Version Control (Git)
- [ ] Complete `assessments/git_assessment.md`
- [ ] Minimum passing: 60% (6/10 questions correct)
- [ ] If you score < 60%: Review `resources/git_primer.md`

---

## Readiness Checklist

After completing all assessments, you should be able to:

### Programming (R or Python)
- [ ] Read and write basic scripts
- [ ] Install packages/libraries
- [ ] Load and manipulate data frames
- [ ] Create simple plots
- [ ] Write functions
- [ ] Use loops and conditionals
- [ ] Debug simple errors

### Command Line
- [ ] Navigate directory structure (cd, ls, pwd)
- [ ] Create, move, copy files and directories
- [ ] View file contents (cat, head, tail, less)
- [ ] Use pipes and redirects (|, >, >>)
- [ ] Run programs with arguments
- [ ] Understand file permissions
- [ ] Use grep for pattern matching

### Biology
- [ ] Explain what RNA is and its role in cells
- [ ] Understand gene expression concept
- [ ] Know difference between DNA, RNA, and protein
- [ ] Explain cell types and why they differ
- [ ] Understand basics of transcription
- [ ] Know what a genome and transcriptome are

### Statistics
- [ ] Understand mean, median, standard deviation
- [ ] Know what a p-value represents
- [ ] Understand correlation vs causation
- [ ] Interpret distributions and histograms
- [ ] Understand concept of statistical significance
- [ ] Know basics of hypothesis testing

### Version Control
- [ ] Clone a repository
- [ ] Pull latest changes
- [ ] Commit changes with meaningful messages
- [ ] Push changes to remote
- [ ] Basic understanding of branches

---

## Recommended Preparation Paths

### If You're Coming From...

#### **Biology Background, Limited Programming**
**Estimated prep time:** 20-40 hours

Priority order:
1. Programming primer (focus on R or Python)
2. Command line basics
3. Git fundamentals
4. Light statistics review

**Start with:** `resources/programming_primer.md`

#### **Computer Science Background, Limited Biology**
**Estimated prep time:** 10-20 hours

Priority order:
1. Biology fundamentals
2. Statistics review (focus on biological applications)
3. Practice with biological datasets

**Start with:** `resources/biology_stats_primer.md`

#### **Wet Lab Biologist, No Computational Experience**
**Estimated prep time:** 40-60 hours

Priority order:
1. Programming fundamentals (choose R or Python)
2. Command line intensive practice
3. Git basics
4. Statistics refresher

**Start with:** Complete Programming Primer, then Command Line Primer

**Recommendation:** Consider taking an introductory programming course first (see External Resources below)

#### **Statistician, New to Biology**
**Estimated prep time:** 15-25 hours

Priority order:
1. Molecular biology fundamentals
2. Genomics-specific concepts
3. RNA-seq data characteristics
4. Practice with biological interpretation

**Start with:** `resources/biology_stats_primer.md`

---

## Self-Assessment Quiz

Take this quick 10-minute quiz to estimate your readiness:

### Quick Readiness Test (Score yourself)

**Programming (4 points max)**
- Can you write a script to read a CSV file? (1 point)
- Can you filter data based on conditions? (1 point)
- Can you create a scatter plot? (1 point)
- Can you write a simple function? (1 point)

**Command Line (3 points max)**
- Can you navigate directories without GUI? (1 point)
- Can you pipe commands together? (1 point)
- Can you search for patterns in files? (1 point)

**Biology (2 points max)**
- Do you know what RNA is and why cells have it? (1 point)
- Do you understand why different cell types exist? (1 point)

**Statistics (1 point max)**
- Do you understand hypothesis testing basics? (1 point)

**TOTAL: _____ / 10 points**

**Interpretation:**
- **8-10 points:** Ready to start! Review primers for any weak areas.
- **5-7 points:** Need focused study in 1-2 areas. 1-2 weeks prep recommended.
- **2-4 points:** Significant prep needed. 4-6 weeks recommended before starting.
- **0-1 points:** Build foundational skills first. Consider introductory courses.

---

## External Resources (Free)

### Programming
- **Python:** [Codecademy Python Course](https://www.codecademy.com/learn/learn-python-3) (25 hours)
- **R:** [Swirl Interactive R Tutorial](https://swirlstats.com/) (10-15 hours)
- **General:** [DataCamp Free Intro](https://www.datacamp.com/courses/free-introduction-to-r) (4 hours)

### Command Line
- [Linux Command Line Basics](https://ubuntu.com/tutorials/command-line-for-beginners) (2 hours)
- [Command Line Crash Course](https://learnpythonthehardway.org/book/appendixa.html) (3 hours)

### Biology & Genomics
- [Khan Academy: Molecular Biology](https://www.khanacademy.org/science/biology/gene-expression-central-dogma) (5 hours)
- [DNA Learning Center: Gene Expression](https://dnalc.cshl.edu/) (Interactive)
- [NIH: Introduction to Genomics](https://www.genome.gov/About-Genomics/Introduction-to-Genomics) (Reading)

### Statistics
- [Khan Academy: Statistics & Probability](https://www.khanacademy.org/math/statistics-probability) (15-20 hours relevant sections)
- [Seeing Theory: Visual Statistics](https://seeing-theory.brown.edu/) (Interactive, 3-5 hours)

### Git
- [GitHub Learning Lab](https://lab.github.com/) (3-5 hours)
- [Git Immersion](https://gitimmersion.com/) (3 hours)

---

## Setting Up Your Environment

Once you've passed the assessments, set up your computational environment:

### Required Software

#### Option 1: Python-Based Workflow
- [ ] Python 3.10+ ([Anaconda](https://www.anaconda.com/products/distribution) recommended)
- [ ] Jupyter Lab or Jupyter Notebook
- [ ] Key packages: `scanpy`, `anndata`, `pandas`, `numpy`, `matplotlib`, `seaborn`

#### Option 2: R-Based Workflow
- [ ] R 4.0+ ([CRAN](https://cran.r-project.org/))
- [ ] RStudio Desktop ([Download](https://www.rstudio.com/products/rstudio/download/))
- [ ] Key packages: `Seurat`, `SingleCellExperiment`, `ggplot2`, `dplyr`

#### Both Paths Need:
- [ ] Git ([Download](https://git-scm.com/downloads))
- [ ] Text editor (VS Code, Sublime, or built-in with IDE)
- [ ] Terminal/Command line access

### Test Your Setup

Run `assessments/environment_test.md` to verify:
- [ ] All required software installed
- [ ] Packages/libraries load correctly
- [ ] Test datasets can be read
- [ ] Plots can be generated
- [ ] Git is configured

---

## What to Expect from the Full Curriculum

### Time Commitment
- **Minimum:** 90-110 hours over 16-20 weeks
- **Realistic with job/studies:** 4-6 months
- **Comfortable pace:** 6-9 months

### Difficulty Curve
```
Easy → → → → → → → → → → Hard
|         |         |         |
Course 1  Course 3  Course 5  Course 6
          Course 2  Course 4
```

### Computational Requirements
- **Minimum:** 8GB RAM, 50GB free disk, 4 CPU cores
- **Recommended:** 16GB RAM, 100GB free disk, 8 CPU cores
- **Some modules:** Cloud computing alternatives provided

### Support Expectations
- Self-paced learning (no fixed schedule)
- Access to course materials 24/7
- Discussion forum for questions (if available)
- Office hours (if offered)

---

## Ready to Start?

### If You Passed All Assessments (≥70% each):
✅ **Proceed to Course 1:** `../1.single_cell_processing_course/START_HERE.md`

### If You Need More Preparation:
📚 **Focus on your weak areas:**
- Review relevant primers in `resources/`
- Complete recommended external courses
- Retake assessments when ready
- No time limit - go at your own pace!

### Questions or Concerns?
- Review `FAQ.md` for common questions
- Check your institution's support resources
- Consider study groups or peer learning

---

## Encouragement

Remember:
- **Everyone starts somewhere** - Even experts were beginners once
- **Progress over perfection** - You don't need to be perfect to start
- **Learning is iterative** - You'll deepen understanding as you practice
- **Ask for help** - No question is too basic
- **Celebrate small wins** - Each completed module is progress!

**Good luck with your single-cell RNA-seq journey! 🧬🔬📊**
