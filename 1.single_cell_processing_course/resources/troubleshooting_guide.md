# Troubleshooting Guide: Course 1 - Single-Cell Processing

**Common errors, their causes, and solutions**

---

## Quick Reference: Error Types

| Error Category | Symptoms | Jump To |
|----------------|----------|---------|
| Installation | Package won't install, version conflicts | [Section 1](#1-installation-issues) |
| File/Path | "File not found", permission denied | [Section 2](#2-file-and-path-issues) |
| Memory | "Out of memory", process killed | [Section 3](#3-memory-issues) |
| Cell Ranger | Alignment fails, low mapping rate | [Section 4](#4-cell-ranger-issues) |
| Data Loading | Can't read matrix, format errors | [Section 5](#5-data-loading-issues) |
| QC/Filtering | Unexpected results, all cells filtered | [Section 6](#6-qc-and-filtering-issues) |
| Doublet Detection | Tool fails, unrealistic scores | [Section 7](#7-doublet-detection-issues) |

---

## 1. Installation Issues

### Problem 1.1: `scanpy` won't install

**Error message:**
```
ERROR: Could not find a version that satisfies the requirement scanpy
```

**Cause:** Python version too old or pip outdated

**Solution:**
```bash
# Check Python version (need 3.7+)
python --version

# Upgrade pip
pip install --upgrade pip

# Try installing again
pip install scanpy

# If still fails, use conda
conda install -c conda-forge scanpy
```

---

### Problem 1.2: Version conflicts with dependencies

**Error message:**
```
ERROR: pip's dependency resolver does not currently take into account all the packages
 that are installed. This behaviour is the source of the following dependency conflicts.
```

**Cause:** Incompatible package versions

**Solution:**
```bash
# Create a fresh environment
conda create -n scrna python=3.10
conda activate scrna

# Install core packages together
pip install scanpy leidenalg python-igraph

# Or use conda for everything
conda install -c conda-forge scanpy leidenalg
```

---

### Problem 1.3: Cell Ranger won't run

**Error message:**
```
cellranger: command not found
```

**Cause:** Cell Ranger not in PATH

**Solution:**
```bash
# Download Cell Ranger from 10x Genomics website
# Extract and add to PATH
export PATH=/path/to/cellranger-8.0.0:$PATH

# Add to .bashrc or .bash_profile for permanent access
echo 'export PATH=/path/to/cellranger-8.0.0:$PATH' >> ~/.bashrc
source ~/.bashrc

# Verify
cellranger --version
```

---

## 2. File and Path Issues

### Problem 2.1: "File not found" error

**Error message:**
```
FileNotFoundError: [Errno 2] No such file or directory: 'data/matrix.mtx'
```

**Cause:** Incorrect path or file doesn't exist

**Solution:**
```python
import os

# Check if file exists
file_path = 'data/matrix.mtx'
print(os.path.exists(file_path))  # Should print True

# Check current directory
print(os.getcwd())

# List files in directory
print(os.listdir('data/'))

# Use absolute path if needed
file_path = os.path.abspath('data/matrix.mtx')
```

---

### Problem 2.2: Permission denied

**Error message:**
```
PermissionError: [Errno 13] Permission denied: '/output/results.h5ad'
```

**Cause:** No write permission in directory

**Solution:**
```bash
# Check permissions
ls -l /output/

# Change permissions (if you own the directory)
chmod 755 /output/

# Or save to a directory you own
mkdir -p ~/my_results
# Then use ~/my_results/ in your code
```

---

### Problem 2.3: Spaces in file paths

**Error message:**
```
FileNotFoundError: [Errno 2] No such file or directory: '/path/My'
```

**Cause:** Path with spaces not properly quoted

**Solution:**
```python
# DON'T: path with spaces unquoted
path = /path/My Documents/file.txt  # WRONG

# DO: Use quotes
path = "/path/My Documents/file.txt"

# Or avoid spaces entirely (rename)
mv "My Documents" My_Documents
```

---

## 3. Memory Issues

### Problem 3.1: Out of memory during Cell Ranger

**Error message:**
```
[error] Pipestance failed. Error log at: ...
...MemoryError...
```

**Cause:** Insufficient RAM for dataset size

**Solution:**
```bash
# Check available memory
free -h

# Option 1: Reduce threads (uses less parallel memory)
cellranger count --localmem=32 --localcores=4 ...

# Option 2: Use --maxjobs to limit parallelism
cellranger count --maxjobs=2 ...

# Option 3: Run on HPC or cloud with more RAM

# Option 4: Downsample FASTQs first (for testing)
seqtk sample -s100 reads_R1.fastq.gz 0.1 > subset_R1.fastq.gz
```

---

### Problem 3.2: Kernel dies when loading data

**Error message:**
```
Kernel Restarted
The kernel for /path/notebook.ipynb appears to have died.
```

**Cause:** Loading large objects into memory

**Solution:**
```python
# Option 1: Load in backed mode (doesn't load into memory)
import scanpy as sc
adata = sc.read_h5ad('large_file.h5ad', backed='r')

# Option 2: Subsample for exploration
adata = sc.read_h5ad('large_file.h5ad')
adata_sub = sc.pp.subsample(adata, n_obs=10000, copy=True)

# Option 3: Use sparse matrices (should be default, but check)
print(adata.X)  # Should say <sparse matrix>

# Option 4: Close other programs to free memory
```

---

### Problem 3.3: "Killed" message during processing

**Error message:**
```
Killed
```

**Cause:** Operating system killed process due to memory usage

**Solution:**
```python
# Process in chunks
import scanpy as sc

# Instead of loading all at once:
# adata = sc.read_h5ad('huge_file.h5ad')

# Load and process in chunks
for chunk_start in range(0, total_cells, chunk_size):
    chunk_end = min(chunk_start + chunk_size, total_cells)
    adata_chunk = sc.read_h5ad('huge_file.h5ad')[chunk_start:chunk_end]
    # Process chunk...
    
# Or use cloud computing (Google Colab, AWS, etc.)
```

---

## 4. Cell Ranger Issues

### Problem 4.1: Very low mapping rate (< 50%)

**Error message:** (in metrics_summary.csv)
```
Reads Mapped to Genome: 25%
```

**Cause:** Wrong reference genome or contamination

**Solution:**
```bash
# Check reference genome matches species
# Human: GRCh38
# Mouse: mm10 or mm39

# Verify FASTQ quality
fastqc sample_S1_L001_R1_001.fastq.gz

# Check for contamination (blast a few reads)

# Try different reference or rebuild
cellranger mkref --genome=GRCh38 \
                 --fasta=GRCh38.fa \
                 --genes=genes.gtf
```

---

### Problem 4.2: No cells detected

**Error message:** (in web_summary.html)
```
Estimated Number of Cells: 0
```

**Cause:** Wrong chemistry, low quality library, or incorrect FASTQ structure

**Solution:**
```bash
# Check chemistry
cellranger count --chemistry=SC3Pv3 ...  # Try specifying explicitly

# Or let Cell Ranger auto-detect
cellranger count --chemistry=auto ...

# Check FASTQ structure
zcat sample_R1.fastq.gz | head -n 20

# Ensure correct FASTQ naming
# Should be: [Sample Name]_S1_L001_R[1/2]_001.fastq.gz

# Check if reads are labeled correctly (index vs R1 vs R2)
```

---

### Problem 4.3: Cell Ranger stuck/hanging

**Error message:** None - process just hangs

**Cause:** File system issues or resource contention

**Solution:**
```bash
# Check if it's really stuck (might just be slow)
tail -f cellranger_output/_log

# Check system resources
top
df -h  # Check disk space

# Kill and restart
pkill -9 cellranger

# Clean up previous partial run
rm -rf cellranger_output/

# Restart with fresh run
cellranger count ...
```

---

## 5. Data Loading Issues

### Problem 5.1: Can't read 10x MTX files

**Error message:**
```
ValueError: cannot open file at 'filtered_feature_bc_matrix/matrix.mtx.gz'
```

**Cause:** Incorrect directory structure or missing files

**Solution:**
```python
import scanpy as sc
import os

# Check directory structure
# Should have: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
path = 'filtered_feature_bc_matrix/'
print(os.listdir(path))

# Correct way to load
adata = sc.read_10x_mtx(
    path,
    var_names='gene_symbols',  # or 'gene_ids'
    cache=True
)

# If files are not gzipped (older format)
# Might need to specify explicitly
```

---

### Problem 5.2: H5AD file won't load

**Error message:**
```
KeyError: 'Unable to open object (object 'X' doesn't exist)'
```

**Cause:** Corrupted file or version incompatibility

**Solution:**
```python
# Check anndata version
import anndata
print(anndata.__version__)

# Try updating anndata
# pip install --upgrade anndata

# If file is corrupted, re-export from original data
# If you have the original AnnData object:
adata.write('file_repaired.h5ad', compression='gzip')

# Check file integrity
import h5py
with h5py.File('file.h5ad', 'r') as f:
    print(list(f.keys()))
```

---

### Problem 5.3: Gene names vs Gene IDs confusion

**Error message:**
```
KeyError: 'TP53' (when trying to plot a gene)
```

**Cause:** Using gene symbols but data has ENSEMBL IDs

**Solution:**
```python
import scanpy as sc

# Check what's in var names
print(adata.var_names[:10])
# If you see ENSG... instead of TP53, etc.

# Option 1: Make gene symbols the index
adata.var_names_make_unique()  # In case of duplicates
adata.var['gene_ids'] = adata.var_names  # Save IDs
adata.var_names = adata.var['gene_symbols']  # Use symbols
adata.var_names_make_unique()

# Option 2: Use the ID directly
# Find the ENSEMBL ID for your gene
tp53_id = adata.var[adata.var['gene_symbols'] == 'TP53'].index[0]
sc.pl.umap(adata, color=tp53_id)
```

---

## 6. QC and Filtering Issues

### Problem 6.1: All or most cells filtered out

**Error message:**
```
After filtering: 234 cells remaining (started with 5000)
```

**Cause:** Too strict QC thresholds

**Solution:**
```python
import scanpy as sc
import numpy as np

# Check QC distributions BEFORE filtering
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])

# Use adaptive thresholds (MAD method)
def is_outlier(adata, metric, nmads):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * np.median(np.abs(M - np.median(M)))) | \
              (np.median(M) + nmads * np.median(np.abs(M - np.median(M))) < M)
    return outlier

# Use gentler thresholds
adata = adata[~is_outlier(adata, 'log1p_total_counts', 5)]  # 5 MADs instead of 3

# Or use tissue-specific fixed thresholds
# Brain: max_mt = 20
# Blood: max_mt = 10
```

---

### Problem 6.2: High mt% cells but they're real

**Error message:** None, but removing high mt% removes expected cell types

**Cause:** Tissue has naturally high mt% (e.g., neurons)

**Solution:**
```python
# Don't use mt% cutoff blindly
# Check if high-mt cells have specific markers

# Subset high mt cells
high_mt = adata[adata.obs['pct_counts_mt'] > 20].copy()

# Check for cell type markers
sc.tl.rank_genes_groups(high_mt, groupby='leiden')
sc.pl.rank_genes_groups(high_mt)

# If they're real cells (have markers), keep them!
# Use higher mt threshold for this tissue
max_mt_threshold = 25  # Instead of 10-15
```

---

### Problem 6.3: QC metrics drive clustering

**Error message:** None, but clusters separate by QC not biology

**Cause:** Under-filtering or batch effects

**Solution:**
```python
# Check if clustering is QC-driven
sc.pl.umap(adata, color=['leiden', 'total_counts', 'pct_counts_mt'])

# If clusters follow QC metrics:
# Option 1: Filter more strictly
adata = adata[adata.obs['pct_counts_mt'] < 15]

# Option 2: Regress out QC effects (use with caution!)
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# Option 3: Re-check for low-quality cells/doublets
# Run Scrublet again with adjusted parameters
```

---

## 7. Doublet Detection Issues

### Problem 7.1: Scrublet predicts all cells as doublets

**Error message:**
```
Detected doublet rate: 95%
```

**Cause:** Wrong expected doublet rate parameter

**Solution:**
```python
import scrublet as scr

# Default expected rate might be too high
scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)  # Try lower rate

# Or use automatic detection
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2,
    min_cells=3,
    min_gene_variability_pctl=85,
    n_prin_comps=30
)

# Check threshold
scrub.plot_histogram()  # Look at distribution
```

---

### Problem 7.2: Scrublet predicts no doublets

**Error message:**
```
Detected doublet rate: 0%
```

**Cause:** Threshold too high or homogeneous cell types

**Solution:**
```python
# Lower threshold
scrub.call_doublets(threshold=0.25)  # Instead of default

# Check simulation vs observed
scrub.plot_histogram()

# If all cell types are similar, doublet detection is hard
# Manually check for high UMI/gene count cells
high_counts = adata[adata.obs['total_counts'] > np.percentile(adata.obs['total_counts'], 95)]
sc.tl.leiden(high_counts, resolution=0.5)
# Check if these cluster by mixed markers
```

---

### Problem 7.3: Doublet detection fails with error

**Error message:**
```
ValueError: could not convert string to float: 'NA'
```

**Cause:** NaN values in count matrix

**Solution:**
```python
# Check for NaN
import numpy as np
print(np.any(np.isnan(adata.X.data)))

# Remove cells/genes with NaN
adata = adata[~np.isnan(adata.X.sum(axis=1))]

# Or fill NaN with 0
adata.X.data[np.isnan(adata.X.data)] = 0

# Ensure sparse matrix
from scipy.sparse import issparse, csr_matrix
if not issparse(adata.X):
    adata.X = csr_matrix(adata.X)
```

---

## 8. General Python/R Errors

### Problem 8.1: ModuleNotFoundError

**Error message:**
```
ModuleNotFoundError: No module named 'scanpy'
```

**Cause:** Package not installed in current environment

**Solution:**
```bash
# Check which Python is being used
which python

# Check environment
conda env list

# Activate correct environment
conda activate scrna

# Install package in current environment
pip install scanpy

# Or check in Jupyter notebook
import sys
print(sys.executable)
```

---

### Problem 8.2: Jupyter kernel issues

**Error message:**
```
Kernel dead or not found
```

**Cause:** Kernel not installed for environment

**Solution:**
```bash
# Install ipykernel in your environment
conda activate scrna
pip install ipykernel

# Register kernel with Jupyter
python -m ipykernel install --user --name=scrna --display-name="scRNA-seq"

# Restart Jupyter and select the correct kernel
```

---

### Problem 8.3: Plots not showing

**Error message:** None, but no plots appear

**Cause:** Wrong backend or inline magic missing

**Solution:**
```python
# In Jupyter notebook, add at top:
%matplotlib inline

# Or for interactive plots:
%matplotlib widget

# Check plotting backend
import matplotlib
print(matplotlib.get_backend())

# Force inline backend
matplotlib.use('inline')

# For scanpy, set DPI
import scanpy as sc
sc.set_figure_params(dpi=80, dpi_save=300)
```

---

## 9. Performance Issues

### Problem 9.1: Code running very slowly

**Symptoms:** Processing takes hours instead of minutes

**Causes & Solutions:**
```python
# Cause 1: Not using sparse matrices
print(type(adata.X))  # Should be sparse
from scipy.sparse import csr_matrix
adata.X = csr_matrix(adata.X)

# Cause 2: Too many threads
# Some operations are slower with many cores
import os
os.environ['OPENBLAS_NUM_THREADS'] = '4'
os.environ['MKL_NUM_THREADS'] = '4'

# Cause 3: Not using optimized libraries
# Install intel math kernel library
conda install mkl

# Cause 4: Recalculating instead of storing
# Store intermediate results
adata.write('checkpoint.h5ad')
```

---

## 10. Reproducibility Issues

### Problem 10.1: Results change every run

**Symptoms:** Different UMAP layout, different clusters

**Cause:** Random seed not set

**Solution:**
```python
# Set random seed
import numpy as np
import random

np.random.seed(42)
random.seed(42)

# For scanpy specifically
import scanpy as sc
sc.settings.seed = 42

# Some algorithms still have randomness; check parameters
sc.tl.leiden(adata, random_state=42)
sc.tl.umap(adata, random_state=42)
```

---

## Need More Help?

### Before Asking for Help:

1. **Check Error Message Carefully**
   - Read the full error, not just the last line
   - Note the line number where error occurred

2. **Search Online**
   - Copy exact error message to Google
   - Check GitHub issues for the package
   - Search Bioconductor/Biostars forums

3. **Check Software Versions**
   ```python
   import scanpy as sc
   import anndata
   print(f"scanpy=={sc.__version__}")
   print(f"anndata=={anndata.__version__}")
   
   # Full environment
   sc.logging.print_versions()
   ```

4. **Create Minimal Reproducible Example**
   - Isolate the problem
   - Use small test data if possible
   - Remove unnecessary code

### Where to Ask:

- **Course Discussion Forum:** [Link]
- **Scanpy GitHub Issues:** https://github.com/scverse/scanpy/issues
- **Bioconductor Support:** https://support.bioconductor.org/
- **Stack Overflow:** Tag with `scanpy`, `single-cell`, etc.

### When Posting:

Include:
- Exact error message (full traceback)
- Software versions (`sc.logging.print_versions()`)
- Minimal code to reproduce
- What you've already tried

---

**Last Updated:** February 2, 2026  
**Contribute:** Found a solution to an unlisted problem? Add it to the course repository!
