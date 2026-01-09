# Environment Setup

## Required Software

### Option 1: Conda Environment (Recommended)

```bash
# Create environment
conda create -n scrna-course python=3.11 -y
conda activate scrna-course

# Core Python packages
pip install scanpy anndata pandas numpy matplotlib seaborn
pip install scrublet

# Jupyter for labs
pip install jupyterlab
```

### Option 2: Docker

```bash
docker pull scanpy/scanpy:latest
```

---

## Bioinformatics Tools

### Cell Ranger (10x Genomics)

1. Download from: https://www.10xgenomics.com/support/software/cell-ranger/downloads
2. Requires registration
3. Extract and add to PATH:

```bash
export PATH=/path/to/cellranger-8.0.0:$PATH
```

### kallisto + bustools

```bash
# Conda install
conda install -c bioconda kallisto bustools

# Or from source
# https://www.kallistobus.tools/
```

### STARsolo

```bash
conda install -c bioconda star
```

---

## Demo Data Download

### 1k PBMCs (Primary dataset for labs)

```bash
mkdir -p data/1k_pbmc
cd data/1k_pbmc

# Download from 10x Genomics
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar

# Download pre-built reference
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```

### 5k PBMCs (Final project option)

```bash
mkdir -p data/5k_pbmc
cd data/5k_pbmc

wget https://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar
tar -xvf 5k_pbmc_v3_fastqs.tar
```

---

## Verify Installation

```bash
# Check Cell Ranger
cellranger --version

# Check kallisto
kallisto version

# Check Python environment
python -c "import scanpy; print(scanpy.__version__)"
```

---

## Recommended Directory Structure

```
single_cell_processing_course/
├── data/
│   ├── 1k_pbmc/
│   │   ├── fastqs/
│   │   └── reference/
│   └── 5k_pbmc/
├── labs/
├── assignments/
├── outputs/
│   ├── cellranger/
│   └── kallisto/
└── resources/
```

Create this structure:

```bash
mkdir -p data/{1k_pbmc,5k_pbmc} outputs/{cellranger,kallisto}
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Cell Ranger memory error | Need ≥16GB RAM, use `--localmem` flag |
| kallisto index fails | Check reference FASTA integrity |
| Python import errors | Verify conda environment is active |

