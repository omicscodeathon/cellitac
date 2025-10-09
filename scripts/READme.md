# scATACtf Setup Guide
## this is th first step before running both R or Python Scripts..

## System Requirements

### Python Environment
**Python Version:** 3.10.17

### R Environment
**R Version:** 4.5.1 or higher

---

## Installation Instructions

### 1. Python Setup

#### Step 1: Create Virtual Environment (Recommended)
```bash
# Create virtual environment
python -m venv scatactf_env

# Activate virtual environment
# On Windows:
scatactf_env\Scripts\activate
# On macOS/Linux:
source scatactf_env/bin/activate
```

#### Step 2: Install Python Dependencies
```bash
pip install -r requirements.txt
```

#### Manual Installation (Alternative)
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

### 2. R Setup

#### Step 1: Install R Packages
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

#### Alternative: Run Installation Script
```r
# From R console
source("install_R_packages.R")
```

---

## Python Package Versions

| Package | Version |
|---------|---------|
| pandas | 2.2.3 |
| numpy | 2.1.2 |
| scikit-learn | 1.6.1 |
| xgboost | 3.0.5 |
| imbalanced-learn | 0.14.0 |
| matplotlib | 3.10.1 |
| seaborn | 0.13.2 |
| plotly | 6.3.0 |
| networkx | 3.3 |
| openpyxl | 3.1.5 |
| joblib | 1.4.2 |
| jupyter | 1.1.1 |

---

## R Package Versions

| Package | Version |
|---------|---------|
| Seurat | 5.3.0 |
| Signac | 1.15.0 |

---

## Verification

### Verify Python Installation
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

print("All Python packages imported successfully!")
```

### Verify R Installation
```r
library(Seurat)
library(Signac)

cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Signac version:", as.character(packageVersion("Signac")), "\n")
```

---

## Troubleshooting

### Python Issues
- **Issue:** pip install fails
  - **Solution:** Try upgrading pip: `pip install --upgrade pip`
  
- **Issue:** Version conflicts
  - **Solution:** Use a fresh virtual environment

### R Issues
- **Issue:** Package installation fails
  - **Solution:** Update R to version 4.5.1 or higher
  
- **Issue:** Bioconductor packages needed
  - **Solution:** Install BiocManager first:
    ```r
    install.packages("BiocManager")
    ```

---

## Contact & Support

For issues or questions, please refer to the project repository or contact the development team.

---

**Note:** The complete setup files (requirements.txt and install_R_packages.R) can be found in the project repository's scripts directory.
