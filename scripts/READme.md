## System Requirements

- **OS:** Linux or macOS (Windows not supported)
- **Python:** 3.9 – 3.12 (not 3.13+)
- **R:** 4.4.3 or higher
- **Conda:** required for bioconda installation (Miniconda or Anaconda)

---

## Installation

### Option 1: bioconda (recommended — installs Python + R dependencies automatically)

```bash
conda install -c bioconda -c conda-forge cellitac
```

### Option 2: PyPI (Python only — R packages must be installed separately)

```bash
pip install cellitac
```

If using PyPI, install R packages manually:

```r
install.packages("Seurat")
install.packages("Signac")
install.packages("hdf5r")
install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("EnsDb.Hsapiens.v75")
```

---

## Quick Start

```bash
cellitac --help
cellitac-preprocess --help
cellitac-model --help
```
