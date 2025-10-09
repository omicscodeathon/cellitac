# scATACtf Setup Guide

> ** IMPORTANT:** Complete this setup process before running any R or Python scripts.

---

## Table of Contents
- [System Requirements](#system-requirements)
- [Installation Instructions](#installation-instructions)
  - [Step 1: R Environment Setup](#step-1-r-environment-setup)
  - [Step 2: Python Environment Setup](#step-2-python-environment-setup)
- [Verification](#verification)
- [Troubleshooting](#troubleshooting)
- [Next Steps](#next-steps)

---

## System Requirements

### R Environment
**Required Version:** R 4.5.1 or higher

### Python Environment
**Required Version:** Python 3.10.17

---

## Installation Instructions

### Step 1: R Environment Setup

#### Install R Packages

Open R or RStudio and run:

```r
# Install required packages
install.packages("Seurat")
install.packages("Signac")

# Verify installations
library(Seurat)
library(Signac)

# Check versions
packageVersion("Seurat")   # Should be 5.3.0
packageVersion("Signac")   # Should be 1.15.0
```

#### Alternative: Use Installation Script

```r
# From R console
source("install_R_packages.R")
```

#### Required R Packages

| Package | Version |
|---------|---------|
| Seurat | 5.3.0 |
| Signac | 1.15.0 |

---

### Step 2: Python Environment Setup

#### Option A: Quick Installation (Recommended)

**1. Create Virtual Environment:**
```bash
# Create virtual environment
python -m venv scatactf_env

# Activate virtual environment
# On Windows:
scatactf_env\Scripts\activate

# On macOS/Linux:
source scatactf_env/bin/activate
```

**2. Install All Dependencies:**
```bash
pip install -r Python_requirements.txt
```

---

#### Option B: Manual Installation

If you prefer to install packages individually:

```bash
pip install pandas==2.2.3
pip install numpy==2.1.2
pip install scikit-learn==1.6.1
pip install xgboost==3.0.5
pip install imbalanced-learn==0.14.0
pip install matplotlib==3.10.1
pip install seaborn==0.13.2
pip install plotly==6.3.0
pip install networkx==3.3
pip install openpyxl==3.1.5
pip install joblib==1.4.2
pip install jupyter==1.1.1
```

---

#### Required Python Packages

| Package | Version | Purpose |
|---------|---------|---------|
| pandas | 2.2.3 | Data manipulation |
| numpy | 2.1.2 | Numerical computing |
| scikit-learn | 1.6.1 | Machine learning |
| xgboost | 3.0.5 | Gradient boosting |
| imbalanced-learn | 0.14.0 | Class imbalance handling |
| matplotlib | 3.10.1 | Visualization |
| seaborn | 0.13.2 | Statistical visualization |
| plotly | 6.3.0 | Interactive plots |
| networkx | 3.3 | Network analysis |
| openpyxl | 3.1.5 | Excel file handling |
| joblib | 1.4.2 | Model serialization |
| jupyter | 1.1.1 | Interactive notebooks |

---

## Verification

### Verify R Installation

Run the following in R/RStudio:

```r
library(Seurat)
library(Signac)

cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Signac version:", as.character(packageVersion("Signac")), "\n")
```

**Expected Output:**
```
Seurat version: 5.3.0
Signac version: 1.15.0
```

---

### Verify Python Installation

Run the following in Python:

```python
import pandas as pd
import numpy as np
import sklearn
import xgboost as xgb
from imblearn.over_sampling import BorderlineSMOTE
import matplotlib.pyplot as plt
import seaborn as sns
import plotly
import networkx as nx

print("✓ All Python packages imported successfully!")
```

**Expected Output:**
```
✓ All Python packages imported successfully!
```

---

## Troubleshooting

### Python Issues

#### Issue: `pip install` fails
**Solution:** Try upgrading pip first
```bash
pip install --upgrade pip
```

#### Issue: Version conflicts
**Solution:** Use a fresh virtual environment
```bash
# Deactivate current environment
deactivate

# Remove old environment
rm -rf scatactf_env  # Linux/macOS
# OR
rmdir /s scatactf_env  # Windows

# Create fresh environment
python -m venv scatactf_env
```

#### Issue: Package installation takes too long
**Solution:** Use a different pip mirror or install specific packages first
```bash
pip install --upgrade pip setuptools wheel
pip install -r Python_requirements.txt
```

---

### R Issues

#### Issue: Package installation fails
**Solution:** Update R to version 4.5.1 or higher
- Download from: https://cran.r-project.org/

#### Issue: Bioconductor packages needed
**Solution:** Install BiocManager first
```r
install.packages("BiocManager")
BiocManager::install()
```

#### Issue: Seurat or Signac installation fails
**Solution:** Install dependencies first
```r
install.packages(c("Matrix", "Rcpp", "RcppArmadillo"))
install.packages("Seurat")
install.packages("Signac")
```

---

## Setup Files Location

All setup files can be found in the project repository's `scripts` directory:

- **Python_requirements.txt** - Python package dependencies
- **install_R_packages.R** - R package installation script

---

## Next Steps

After completing the setup:

1.  Verify both R and Python installations (see [Verification](#verification) section)
2.  Proceed to the pipeline workflow guide
3.  Start with data preprocessing (R scripts)
4.  Continue with machine learning analysis (Python scripts)

** For detailed pipeline instructions, see:** [Pipeline Workflow Guide](https://github.com/omicscodeathon/scatactf/blob/main/workflow/pipeline.md)

---

## Contact & Support

For issues or questions, please:
- Open an issue on the GitHub repository
- Contact the development team
- Refer to the documentation

---

**Note:** This setup guide is the **first step** before running any scATACtf analysis scripts. Make sure all packages are installed and verified before proceeding with the pipeline.
