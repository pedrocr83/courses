# Environment Setup Verification

**Time Required:** 15-30 minutes  
**Purpose:** Verify your computational environment is ready for the single-cell RNA-seq curriculum

---

## Instructions

1. Follow each section below in order
2. Run the provided test commands in your terminal
3. Check off each item as it passes
4. If any test fails, see the Troubleshooting section at the bottom

---

## 1. Python Environment

*Skip this section if you are using R only.*

### Python Installation
- [ ] Python 3.10+ installed

```bash
python3 --version
# Expected: Python 3.10.x or higher
```

### Package Manager
- [ ] pip or conda available

```bash
pip --version
# OR
conda --version
```

### Jupyter
- [ ] Jupyter Lab or Notebook installed

```bash
jupyter --version
```

### Key Packages
- [ ] All required Python packages installed

Run this test script:
```python
import sys
print(f"Python: {sys.version}")

packages = {
    'scanpy': 'sc',
    'anndata': 'ad',
    'pandas': 'pd',
    'numpy': 'np',
    'matplotlib': 'mpl',
    'seaborn': 'sns'
}

for pkg_name, alias in packages.items():
    try:
        mod = __import__(pkg_name)
        print(f"  {pkg_name}: {mod.__version__}")
    except ImportError:
        print(f"  {pkg_name}: NOT INSTALLED - run: pip install {pkg_name}")
```

---

## 2. R Environment

*Skip this section if you are using Python only.*

### R Installation
- [ ] R 4.0+ installed

```bash
R --version | head -1
# Expected: R version 4.x.x
```

### RStudio
- [ ] RStudio Desktop installed
- Test: Open RStudio from your applications

### Key Packages
- [ ] All required R packages installed

Run this test script in R:
```r
cat("R version:", R.version.string, "\n\n")

packages <- c("Seurat", "SingleCellExperiment", "ggplot2", "dplyr")

for (pkg in packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, ": ", as.character(packageVersion(pkg)), "\n")
  } else {
    cat(pkg, ": NOT INSTALLED\n")
    cat("  Install with: install.packages('", pkg, "')\n", sep = "")
  }
}
```

---

## 3. Command Line Tools

- [ ] Terminal access works (you're reading this, so probably yes)

### Common Utilities
- [ ] Core tools available

```bash
# Test each tool (should print version or usage info)
grep --version | head -1
awk --version | head -1
sort --version | head -1
```

### Download Tools
- [ ] wget or curl available

```bash
curl --version | head -1
# OR
wget --version | head -1
```

### Bioinformatics Tools (Optional)
- [ ] samtools installed (recommended but not required for Course 1)

```bash
samtools --version | head -1
```

---

## 4. Git

- [ ] Git installed

```bash
git --version
# Expected: git version 2.x.x
```

- [ ] Git configured with name and email

```bash
git config user.name
git config user.email
# Both should print your name and email
```

If not configured:
```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

- [ ] Can clone a repository

```bash
git clone --depth 1 https://github.com/octocat/Hello-World.git /tmp/git-test
rm -rf /tmp/git-test
echo "Git clone: OK"
```

---

## 5. System Resources

### Disk Space
- [ ] At least 50 GB free disk space

```bash
# Linux / macOS
df -h ~ | tail -1

# Windows (PowerShell)
# Get-PSDrive C | Select-Object Free
```

### RAM
- [ ] At least 8 GB RAM (16 GB recommended)

```bash
# Linux
free -h | grep Mem

# macOS
sysctl hw.memsize | awk '{print $2/1024/1024/1024 " GB"}'
```

---

## 6. Quick Functional Test

### Python Functional Test

Save and run this script as `test_environment.py`:

```python
#!/usr/bin/env python3
"""Quick environment test for single-cell RNA-seq curriculum."""

import sys

def test_environment():
    errors = []

    # Test imports
    try:
        import pandas as pd
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
    except ImportError as e:
        errors.append(f"Import error: {e}")
        return errors

    # Test data operations
    try:
        df = pd.DataFrame({
            'gene': ['TP53', 'BRCA1', 'MYC', 'GAPDH', 'ACTB'],
            'expression': [100, 50, 75, 200, 180],
            'pvalue': [0.001, 0.05, 0.01, 0.8, 0.9]
        })
        filtered = df[df['pvalue'] < 0.05]
        assert len(filtered) == 2, f"Expected 2 rows, got {len(filtered)}"
    except Exception as e:
        errors.append(f"Data operation error: {e}")

    # Test plotting
    try:
        fig, ax = plt.subplots()
        ax.scatter(df['expression'], -np.log10(df['pvalue']))
        ax.set_xlabel('Expression')
        ax.set_ylabel('-log10(p-value)')
        ax.set_title('Test Plot')
        fig.savefig('/tmp/env_test_plot.png')
        plt.close()
    except Exception as e:
        errors.append(f"Plotting error: {e}")

    return errors

if __name__ == '__main__':
    errors = test_environment()
    if errors:
        print("ISSUES FOUND:")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print("Environment test PASSED!")
        sys.exit(0)
```

### R Functional Test

Save and run this script as `test_environment.R`:

```r
#!/usr/bin/env Rscript
# Quick environment test for single-cell RNA-seq curriculum

errors <- c()

# Test data operations
tryCatch({
  df <- data.frame(
    gene = c("TP53", "BRCA1", "MYC", "GAPDH", "ACTB"),
    expression = c(100, 50, 75, 200, 180),
    pvalue = c(0.001, 0.05, 0.01, 0.8, 0.9)
  )
  filtered <- df[df$pvalue < 0.05, ]
  stopifnot(nrow(filtered) == 2)
}, error = function(e) {
  errors <<- c(errors, paste("Data operation error:", e$message))
})

# Test plotting
tryCatch({
  png("/tmp/env_test_plot.png")
  plot(df$expression, -log10(df$pvalue),
       xlab = "Expression", ylab = "-log10(p-value)",
       main = "Test Plot", pch = 19)
  dev.off()
}, error = function(e) {
  errors <<- c(errors, paste("Plotting error:", e$message))
})

# Test ggplot2
tryCatch({
  library(ggplot2)
  p <- ggplot(df, aes(x = expression, y = -log10(pvalue))) +
    geom_point() +
    theme_minimal()
  ggsave("/tmp/env_test_ggplot.png", p, width = 5, height = 4)
}, error = function(e) {
  errors <<- c(errors, paste("ggplot2 error:", e$message))
})

if (length(errors) > 0) {
  cat("ISSUES FOUND:\n")
  for (e in errors) cat("  -", e, "\n")
  quit(status = 1)
} else {
  cat("Environment test PASSED!\n")
  quit(status = 0)
}
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `python3: command not found` | Install Python from [python.org](https://www.python.org/) or use Anaconda |
| `pip: command not found` | Run `python3 -m ensurepip --upgrade` |
| Package not found (Python) | `pip install package_name` or `conda install package_name` |
| Package not found (R) | `install.packages("package_name")` in R console |
| Bioconductor packages (R) | `BiocManager::install("package_name")` |
| `git: command not found` | Install from [git-scm.com](https://git-scm.com/) |
| Permission denied | Check file permissions; use `chmod +x script.py` |
| Not enough disk space | Free up space or use an external drive |
| Plotting fails (no display) | Set `matplotlib.use('Agg')` (Python) or use `png()`/`pdf()` device (R) |
| scanpy install fails | Try `pip install scanpy` or `conda install -c conda-forge scanpy` |

### Installation Guides
- **Anaconda (Python):** https://docs.anaconda.com/anaconda/install/
- **R & RStudio:** https://posit.co/download/rstudio-desktop/
- **Git:** https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
- **Bioconductor:** https://www.bioconductor.org/install/

---

## Next Steps

### If All Checks Pass:
✅ Your environment is ready! Proceed to **Course 1:** `../../1.single_cell_processing_course/START_HERE.md`

### If Some Checks Fail:
📚 Follow the troubleshooting steps above, then re-run the failed tests.  
If you're still stuck, consult the installation guides or ask for help.

---

## Checklist Summary

Copy and fill this in:

```
Python Environment:  [ ] Pass  [ ] Skip  [ ] Fail
R Environment:       [ ] Pass  [ ] Skip  [ ] Fail
Command Line Tools:  [ ] Pass  [ ] Fail
Git:                 [ ] Pass  [ ] Fail
System Resources:    [ ] Pass  [ ] Fail
Functional Test:     [ ] Pass  [ ] Fail
```

**Overall Status:** _______________
