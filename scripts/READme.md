# cellitac

[![PyPI version](https://badge.fury.io/py/cellitac.svg)](https://pypi.org/project/cellitac/)
[![Conda Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://anaconda.org/channels/bioconda/packages/cellitac/overview)

**Cell type Identification using Transcription factor Analysis and Chromatin accessibility**

cellitac is a machine learning pipeline for classifying cell types from Single-Cell ATAC + RNA Multiome data using transcription factor analysis and chromatin accessibility features.

> ⚠️ **Before running cellitac, please read the full system requirements and installation instructions below.**

---

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

### Option 2: PyPI (Python only — R packages must be installed manually)

```bash
pip install cellitac
```

If using PyPI, install R packages manually in R:

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

## Python Dependencies

| Package | Version |
|---------|---------|
| Python | 3.11.13 |
| numpy | 2.4.2 |
| pandas | 2.3.3 |
| scikit-learn | 1.8.0 |
| xgboost | 3.2.0 |
| imbalanced-learn | 0.14.1 |
| matplotlib | 3.9.1 |
| seaborn | 0.13.2 |
| plotly | 6.5.2 |
| networkx | 3.6.1 |
| openpyxl | 3.1.5 |
| rpy2 | 3.5.11 |

---

## R Dependencies

| Package | Source |
|---------|--------|
| R | 4.4.3 |
| Seurat | CRAN |
| Signac | CRAN |
| SingleR | Bioconductor |
| celldex | Bioconductor |
| EnsDb.Hsapiens.v75 | Bioconductor |
| hdf5r | CRAN |
| Matrix | CRAN |

---

## Quick Start

```bash
cellitac --help
cellitac-preprocess --help
cellitac-model --help
```

---

## Contributors

- Rana H. Abu-Zeid
- Dr.Olaitan Awe
- Syrus Semawule
- Emmanuel Aroma
- Toheeb Jumah
- Dr.Derek Reiman
