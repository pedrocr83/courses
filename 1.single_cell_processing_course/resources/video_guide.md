# Video Resources Guide

## Recommended Videos by Week & Module

---

## Week 1: Foundations & Raw Data

### Module 1: Introduction to scRNA-seq Processing

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| Single Cell RNA-seq: An Introductory Overview | StatQuest | 10 min | https://www.youtube.com/watch?v=k9VFNLLQP8c |
| What is Single Cell Sequencing? | iBiology | 15 min | https://www.youtube.com/watch?v=3dqXcZBi5vQ |
| 10x Genomics Single Cell Gene Expression | 10x Genomics | 5 min | https://www.youtube.com/watch?v=0Xvxl3i_Ct8 |

**Watch before:** Starting Module 1  
**Key concepts:** Why single-cell matters, platform overview, droplet technology

---

### Module 2: Experimental Design & Metadata

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| Experimental Design for Single-Cell Studies | Broad Institute | 45 min | https://www.youtube.com/watch?v=qDBb0wGMgKk |
| Batch Effects in Single Cell RNA-seq | Harvard Chan Bioinformatics | 20 min | https://www.youtube.com/watch?v=0PmbTlKk3Sg |

**Watch before:** Lab 2  
**Key concepts:** Replicates, batch structure, metadata importance

---

### Module 3: FASTQ Files & Read Structure

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| FASTQ Files Explained | OMGenomics | 8 min | https://www.youtube.com/watch?v=qf4Xv1e-RLk |
| Understanding Phred Quality Scores | Applied Biological Materials | 5 min | https://www.youtube.com/watch?v=1rVDCVZoD9M |
| 10x Chromium Single Cell 3' Gene Expression | 10x Genomics | 8 min | https://www.youtube.com/watch?v=4WiT0W3mj9c |

**Watch before:** Lab 3  
**Key concepts:** FASTQ format, quality scores, 10x read structure

---

## Week 2: Alignment, Quantification & Cell Calling

### Module 4: Reference Preparation

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| Reference Genomes and Gene Annotations | OMGenomics | 12 min | https://www.youtube.com/watch?v=6yDJo4pTQ3Q |
| GTF and GFF Files Explained | Bioinformatics Coach | 10 min | https://www.youtube.com/watch?v=9EkV7FYXbzo |

**Watch before:** Lab 4  
**Key concepts:** FASTA, GTF, genome versions, annotation sources

---

### Module 5: Alignment & Quantification Pipelines

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| Cell Ranger Pipeline Overview | 10x Genomics | 15 min | https://www.youtube.com/watch?v=Z9Yv6R8IgU4 |
| RNA-seq Alignment: STAR Aligner | MIT CompBio | 25 min | https://www.youtube.com/watch?v=hWNWklDGKz8 |
| kallisto and Pseudoalignment | Lior Pachter Lab | 30 min | https://www.youtube.com/watch?v=e4QFXC0B8CI |
| Running Cell Ranger count | 10x Genomics | 12 min | https://www.youtube.com/watch?v=XD0BO1N5wvw |

**Watch before:** Lab 5  
**Key concepts:** Alignment vs pseudoalignment, Cell Ranger workflow, kallisto speed

---

### Module 6: Cell Barcode Detection & UMI Deduplication

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| UMIs in Single Cell RNA-seq | StatQuest | 12 min | https://www.youtube.com/watch?v=xQ6AhGiPjAc |
| EmptyDrops: Distinguishing Cells from Background | Sanger Institute | 20 min | https://www.youtube.com/watch?v=B3oO7l0Lp9k |
| Understanding the Knee Plot | 10x Genomics | 8 min | https://www.youtube.com/watch?v=UVvB6S6Qzz0 |

**Watch before:** Lab 6  
**Key concepts:** UMI purpose, PCR deduplication, cell calling algorithms

---

## Week 3: Quality Control, Filtering & Output

### Module 7: The Count Matrix

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| AnnData: The Single Cell Data Container | Theis Lab | 15 min | https://www.youtube.com/watch?v=jwmuhE9F06A |
| Sparse Matrices Explained | 3Blue1Brown (related) | 10 min | https://www.youtube.com/watch?v=6AoKz_VEUXI |

**Watch before:** Lab 7  
**Key concepts:** AnnData structure, sparse matrices, data organization

---

### Module 8: Quality Control Metrics

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| Single Cell QC: What to Look For | Sanger Institute | 25 min | https://www.youtube.com/watch?v=7t8E6h4kXKE |
| Scanpy Preprocessing Tutorial | Theis Lab | 40 min | https://www.youtube.com/watch?v=uvyG9yLuNSE |
| Mitochondrial Content in scRNA-seq | Broad Institute | 15 min | https://www.youtube.com/watch?v=qPx6y3zPxeA |

**Watch before:** Lab 8  
**Key concepts:** QC metrics, violin plots, threshold selection

---

### Module 9: Filtering Strategies

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| Doublet Detection in Single Cell Data | Satija Lab | 20 min | https://www.youtube.com/watch?v=2dE8P2BsY4A |
| Scrublet Tutorial | swolock | 15 min | https://www.youtube.com/watch?v=0WqH5Kz3FEs |
| QC Filtering Best Practices | EMBL-EBI Training | 30 min | https://www.youtube.com/watch?v=Cx-60EkzZyI |

**Watch before:** Lab 9  
**Key concepts:** Doublet mechanisms, filtering rationale, over-filtering dangers

---

### Module 10: Data Export & Reproducibility

| Video | Source | Duration | Link |
|-------|--------|----------|------|
| Reproducible Bioinformatics with Snakemake | Johannes Köster | 45 min | https://www.youtube.com/watch?v=8xnGFXIc0Uw |
| HDF5 and AnnData File Formats | Scanpy Team | 12 min | https://www.youtube.com/watch?v=3cN-c1RHvKE |

**Watch after:** Completing Lab 10  
**Key concepts:** File formats, reproducibility, pipeline automation

---

## Full Course Playlists

### Comprehensive Single-Cell Analysis

| Playlist | Source | Videos | Link |
|----------|--------|--------|------|
| Single Cell RNA-seq Analysis | Broad Institute | 12 videos | https://www.youtube.com/playlist?list=PLblh5JKOoLUICTaGLRoHQDuF_7q2GfuJF |
| Scanpy Tutorials | Theis Lab | 8 videos | https://www.youtube.com/playlist?list=PLv_9QAaMNdNOqr6Aq2O4h3HxvYYbP3vDN |
| Single Cell Genomics Day | Sanger Institute | 20 videos | https://www.youtube.com/playlist?list=PLiR7IfVl66U_1EJnfQSwJkNEfpULMl0l- |
| 10x Genomics Tutorials | 10x Genomics | 50+ videos | https://www.youtube.com/c/10xGenomics/playlists |

---

## Suggested Viewing Schedule

### Before Course Starts (Pre-work)
- [ ] StatQuest: Single Cell RNA-seq Overview (10 min)
- [ ] iBiology: What is Single Cell Sequencing? (15 min)
- [ ] OMGenomics: FASTQ Files Explained (8 min)

### Week 1
| Day | Video | Duration |
|-----|-------|----------|
| Mon | 10x Genomics Single Cell Overview | 5 min |
| Tue | Experimental Design (Broad) | 45 min |
| Wed | FASTQ + Quality Scores | 13 min |
| Thu | 10x Read Structure | 8 min |
| Fri | *Catch-up / Review* | - |

### Week 2
| Day | Video | Duration |
|-----|-------|----------|
| Mon | Reference Genomes + GTF | 22 min |
| Tue | Cell Ranger Pipeline | 15 min |
| Wed | STAR / kallisto comparison | 55 min |
| Thu | UMIs + Knee Plot | 20 min |
| Fri | EmptyDrops | 20 min |

### Week 3
| Day | Video | Duration |
|-----|-------|----------|
| Mon | AnnData Tutorial | 15 min |
| Tue | Scanpy Preprocessing | 40 min |
| Wed | QC Best Practices | 30 min |
| Thu | Doublet Detection | 35 min |
| Fri | Reproducibility | 45 min |

---

## Video Notes Template

Use this template when watching videos:

```markdown
## Video: [Title]
**Date watched:** 
**Module:** 

### Key Takeaways
1. 
2. 
3. 

### New Terms
- Term: Definition

### Questions to Follow Up
- 

### How This Applies to the Course
- 
```

---

## Additional Resources

### Conference Talks (Advanced)
- Single Cell Genomics Conference recordings
- ISMB/ECCB bioinformatics sessions
- Cold Spring Harbor single-cell courses

### Podcasts
- The Bioinformatics Chat
- Micro Binfie Podcast

### Interactive Tutorials
- Galaxy Training Network: https://training.galaxyproject.org/
- EMBL-EBI Training: https://www.ebi.ac.uk/training/

---

*Note: Video links were current as of course creation. If a link is broken, search for the title on YouTube or the source channel.*

