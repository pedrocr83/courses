# Self-Assessment Checklist: Course 3 - Clustering & Cell Type Annotation

**Track your mastery of dimensionality reduction, clustering, and annotation**

---

## Module 2: Principal Component Analysis

### Concepts
- [ ] I understand why dimensionality reduction is necessary
- [ ] I can explain how PCA works (intuitively)
- [ ] I know what "variance explained" means
- [ ] I understand the difference between analysis PCs and visualization PCs

### Skills
- [ ] I can select HVGs appropriately
- [ ] I can run PCA on scRNA-seq data
- [ ] I can create and interpret an elbow plot
- [ ] I can choose the number of PCs to use

### Can I Explain...
- [ ] Why we use 2000-5000 HVGs?
- [ ] What PC loadings represent?
- [ ] How to read variance explained plots?

**Self-Test:** Generate PCA on PBMC data, justify your choice of nPCs.

---

## Module 3: UMAP & t-SNE

### Concepts
- [ ] I understand non-linear dimensionality reduction
- [ ] I know differences between UMAP and t-SNE
- [ ] I understand key parameters (perplexity, n_neighbors, min_dist)
- [ ] I know UMAP limitations (distances not interpretable)

### Skills
- [ ] I can generate UMAP embeddings
- [ ] I can generate t-SNE embeddings
- [ ] I can tune parameters appropriately
- [ ] I can set random seed for reproducibility

### Red Flags I Can Spot
- [ ] Over-interpreting UMAP distances
- [ ] Using UMAP for quantitative analysis
- [ ] Changing embeddings with minor parameter tweaks

**Self-Test:** Create 3 UMAPs with different parameters, explain differences.

---

## Module 5: Graph-Based Clustering

### Concepts
- [ ] I understand k-NN and SNN graphs
- [ ] I know how Leiden/Louvain algorithms work
- [ ] I understand the resolution parameter
- [ ] I can distinguish over vs under-clustering

### Skills
- [ ] I can build nearest neighbor graphs
- [ ] I can cluster with Leiden algorithm
- [ ] I can test multiple resolutions
- [ ] I can choose appropriate resolution

### Red Flags I Can Spot
- [ ] Too many tiny clusters (over-clustering)
- [ ] Known cell types merged (under-clustering)
- [ ] Clusters with no specific markers

**Self-Test:** Cluster PBMC at resolutions 0.3, 0.8, 1.2—which is best and why?

---

## Module 6: Cluster Evaluation

### Concepts
- [ ] I know cluster quality metrics (silhouette, etc.)
- [ ] I understand biological vs algorithmic validation
- [ ] I know when clusters are "good enough"

### Skills
- [ ] I can calculate silhouette scores
- [ ] I can assess cluster stability
- [ ] I can validate clusters with markers

### Can I Explain...
- [ ] What makes a "good" cluster?
- [ ] When to merge vs split clusters?

**Self-Test:** Evaluate clustering quality, identify problematic clusters.

---

## Module 7: Marker Gene Identification

### Concepts
- [ ] I understand marker gene definition
- [ ] I know different statistical tests for markers (Wilcoxon, t-test, etc.)
- [ ] I understand log2FC and specificity
- [ ] I know how to interpret pct.1 and pct.2

### Skills
- [ ] I can find markers for each cluster
- [ ] I can filter markers by log2FC and FDR
- [ ] I can create marker dotplots
- [ ] I can create marker heatmaps

### Red Flags I Can Spot
- [ ] Markers that aren't specific
- [ ] Ribosomal/mitochondrial "markers"
- [ ] Markers with low effect size

**Self-Test:** Find top 10 markers per cluster, assess their specificity.

---

## Module 8: Manual Annotation

### Concepts
- [ ] I understand annotation confidence levels (HIGH/MEDIUM/LOW)
- [ ] I know canonical markers for major cell types
- [ ] I know where to find marker databases (PanglaoDB, CellMarker)
- [ ] I understand annotation is interpretation, not fact

### Skills
- [ ] I can use marker databases effectively
- [ ] I can assign cell type labels with evidence
- [ ] I can use the Annotation Confidence Framework
- [ ] I can document annotation evidence

### Red Flags I Can Spot
- [ ] Over-confident annotations without evidence
- [ ] Ignoring conflicting markers
- [ ] Lineage-incompatible marker co-expression (doublets)

**Self-Test:** Manually annotate PBMC clusters, assign confidence levels.

---

## Module 9: Automated Annotation

### Concepts
- [ ] I understand reference-based annotation (SingleR)
- [ ] I know marker-based annotation (CellTypist, scType)
- [ ] I understand annotation confidence scores
- [ ] I know when automated methods fail

### Skills
- [ ] I can run SingleR annotation
- [ ] I can run CellTypist (or similar)
- [ ] I can compare manual vs automated
- [ ] I can resolve disagreements

### Red Flags I Can Spot
- [ ] Low confidence scores
- [ ] Disagreement between methods
- [ ] Novel cell types missed by automated tools

**Self-Test:** Run 2 automated tools, compare results, resolve conflicts.

---

## Module 10: Annotation Refinement

### Concepts
- [ ] I understand hierarchical annotation
- [ ] I know when to sub-cluster
- [ ] I know when to merge clusters
- [ ] I understand doublet contamination effects

### Skills
- [ ] I can sub-cluster populations
- [ ] I can merge over-split clusters
- [ ] I can identify doublet clusters
- [ ] I can create final annotation with evidence table

### Red Flags I Can Spot
- [ ] Arbitrary cluster splits
- [ ] Doublet clusters included
- [ ] Inconsistent annotation granularity

**Self-Test:** Identify clusters needing sub-clustering or merging, justify decisions.

---

## Key Competencies Assessment

### Dimensionality Reduction (Modules 2-3)
**Can you:**
- [ ] Select appropriate number of PCs?
- [ ] Generate interpretable UMAP?
- [ ] Explain why we don't use UMAP for quantitative analysis?

**Score:** ___/3

---

### Clustering (Modules 4-6)
**Can you:**
- [ ] Choose appropriate resolution parameter?
- [ ] Identify over/under-clustering?
- [ ] Evaluate cluster quality?
- [ ] Validate clusters with markers?

**Score:** ___/4

---

### Annotation (Modules 7-10)
**Can you:**
- [ ] Find specific marker genes?
- [ ] Assign cell type labels with evidence?
- [ ] Use automated tools appropriately?
- [ ] Assign confidence levels?
- [ ] Document all decisions?

**Score:** ___/5

---

## Overall Mastery Levels

### Basic (Can proceed to Course 4)
- [ ] Can cluster data and get reasonable results
- [ ] Can identify obvious cell types
- [ ] Understands key concepts
- [ ] **Competency scores:** ≥50% each section

### Proficient (Ready for independent analysis)
- [ ] Can optimize clustering systematically
- [ ] Can annotate with appropriate confidence
- [ ] Can troubleshoot problems
- [ ] **Competency scores:** ≥75% each section

### Expert (Can teach/mentor others)
- [ ] Can handle novel cell types
- [ ] Can critically evaluate methods
- [ ] Can integrate multiple lines of evidence
- [ ] **Competency scores:** ≥90% each section

---

## Common Struggles & How to Overcome

### "I can't choose the right resolution"
**Solution:**
1. Test multiple resolutions (0.3, 0.6, 0.9, 1.2)
2. Plot markers at each resolution
3. Choose resolution where known markers separate cleanly
4. Document why you chose it

### "My clusters don't have good markers"
**Red flags this indicates:**
- [ ] Clusters driven by QC not biology → Re-check filtering
- [ ] Over-clustering → Lower resolution
- [ ] Poor data quality → Review Module 8 from Course 1

### "I don't know what this cluster is"
**That's OK!** Options:
- [ ] Label as "Unknown" with LOW confidence
- [ ] Sub-cluster to see if it separates
- [ ] Check for doublets
- [ ] Propose validation experiments

### "Automated tools disagree with manual annotation"
**Investigate:**
- [ ] Check confidence scores (low = unreliable)
- [ ] Verify marker expression manually
- [ ] Consider if automated tool is wrong (they can be!)
- [ ] Report disagreement and your reasoning

---

## Before Assignments

### Assignment 1: Dimensionality Reduction
- [ ] Modules 2-3 mostly checked
- [ ] Can run PCA and UMAP
- [ ] Can compare parameter choices

### Assignment 2: Clustering
- [ ] Modules 4-6 mostly checked
- [ ] Can test multiple resolutions
- [ ] Can evaluate cluster quality

### Assignment 3: Complete Annotation
- [ ] Modules 7-10 mostly checked
- [ ] Can assign confident labels
- [ ] Can document evidence thoroughly

### Final Project
- [ ] ALL modules ≥80% checked
- [ ] Can execute complete workflow
- [ ] Can justify all decisions

---

## Practice Exercises

### Exercise 1: Parameter Exploration
Dataset: PBMC 3k
- [ ] Try n_neighbors = [10, 20, 30, 50]
- [ ] Try resolution = [0.3, 0.6, 0.9, 1.2]
- [ ] Document how results change
- [ ] Choose optimal parameters with justification

### Exercise 2: Marker Validation
- [ ] Find top 5 markers for each cluster
- [ ] Check if they're in PanglaoDB
- [ ] Create dotplot showing expression
- [ ] Identify any questionable markers

### Exercise 3: Annotation Challenge
- [ ] Annotate all clusters manually
- [ ] Run 2 automated tools
- [ ] Compare all 3 annotations
- [ ] Create evidence table
- [ ] Assign confidence levels
- [ ] Write 1-paragraph interpretation

---

## Integration with Other Courses

### From Course 1 (Processing)
**Should be comfortable with:**
- [ ] Loading count matrices
- [ ] QC metrics interpretation
- [ ] Filtering decisions

**If not:** Review Course 1 Modules 7-9

### To Course 4 (Integration)
**Should be ready for:**
- [ ] Multi-batch datasets
- [ ] Batch-aware annotation
- [ ] Integration evaluation

### To Course 2 (DE Analysis)
**Annotations enable:**
- [ ] Cell type-specific DE
- [ ] Pseudobulk aggregation
- [ ] Biological interpretation

---

## Progress Tracking

**Date started:** _______________  
**Current module:** _______________

**Module completion:**
- Module 2 (PCA): ____%
- Module 3 (UMAP/t-SNE): ____%
- Module 5 (Clustering): ____%
- Module 6 (Evaluation): ____%
- Module 7 (Markers): ____%
- Module 8 (Manual annotation): ____%
- Module 9 (Automated): ____%
- Module 10 (Refinement): ____%

**Concepts I need to review:**
1. _______________
2. _______________
3. _______________

---

## Resources for Deeper Learning

### Marker Databases
- [ ] Explored PanglaoDB thoroughly
- [ ] Used CellMarker for my tissue
- [ ] Checked Azimuth references

### Key Papers
- [ ] Read Leiden clustering paper
- [ ] Read UMAP paper (optional but helpful)
- [ ] Read cell type annotation benchmarks

### Practice Datasets
- [ ] PBMC 3k (basic practice)
- [ ] PBMC 10k (more complexity)
- [ ] Tissue-specific dataset relevant to research

---

## Final Readiness Check

Before claiming course completion:

**Technical Skills:**
- [ ] Can independently cluster scRNA-seq data
- [ ] Can optimize parameters systematically
- [ ] Can annotate with appropriate confidence
- [ ] Can create publication-quality figures

**Conceptual Understanding:**
- [ ] Knows when methods fail
- [ ] Can troubleshoot problems
- [ ] Understands limitations
- [ ] Can interpret results biologically

**Documentation:**
- [ ] Creates evidence tables for annotations
- [ ] Justifies all parameter choices
- [ ] Reports confidence levels honestly
- [ ] Produces reproducible analyses

**Ready for Course 4 (Data Integration)?** ✅

---

**Remember:** Good annotation requires patience, critical thinking, and honesty about uncertainty!
