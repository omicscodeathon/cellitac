 # cellitac: Cell type Identification using Transcription factor Analysis and Chromatin accessibility


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![PyPI version](https://badge.fury.io/py/cellitac.svg)](https://pypi.org/project/cellitac/)

<p align="center">
  <img src="https://github.com/omicscodeathon/cellitac/blob/main/figures/logo.png" alt="cellitac logo" width="300"  />
</p>

## Table of Contents
1. [Background](#1-background)
2. [Installation](#2-installation)
3. [Workflow](#3-workflow)
4. [Code Availability](#4-code-avilability)
5. [Reproducibility](#5-reproducibility)
6. [License](#6-license)
7. [Contributors](#7-Contributors)

----
<br>
<p align="center">
  <b>Omics Codeathon General Application - October 2025</b><br>
  Organized by the African Society for Bioinformatics and Computational Biology (ASBCB) with support from the NIH Office of Data Science Strategy.<br>
</p>

----


## 1. Background

Single-cell chromatin accessibility sequencing (scATAC-seq) enables genome-wide profiling of regulatory elements at single-cell resolution. Traditional pipelines identify accessible regions first, then infer TF activity, limiting comprehensive understanding of regulatory programs driving cellular identity. This study develops a robust TF-centric machine learning framework to classify PBMC single-cell datasets using inferred chromVAR transcription factor activities. Our approach addresses data quality challenges through unsupervised redundancy filtering and class-imbalance handling via weighted loss functions, and employs multiple machine learning models for rigorous classification. The resulting computational pipeline enhances single-cell analysis capabilities and provides a systematic approach for discovering TF regulatory networks in immune cell populations.

----

## 2. Installation

The **cellitac** pipeline is designed for easy and direct use through official package managers, ensuring a reproducible environment for single-cell analysis.

### Official Packages
* **PyPI (Python):** [cellitac on PyPI](https://pypi.org/project/cellitac)

---

### 3. Workflow  

<p align="center">
    <img src="https://github.com/omicscodeathon/cellitac/blob/main/figures/cellitac_fixed_methodology.drawio.png" alt="cellitac" width="700" />
</p>  
<p align="center">
    <b>Figure 1.</b> Workflow of the methods employed in this study
</p>



---


## 

### 4. Code Avilability:

All scripts for the **cellitac** project (Python & R) are available in the repository:

👉 Browse the scripts: [Scripts Running](scripts/)

---

###  Demonstration Data  
- Public dataset: [PBMC from a Healthy Donor (10k, 10x Genomics)](https://www.10xgenomics.com/welcome?closeUrl=%2Fdatasets&lastTouchOfferName=PBMC%20from%20a%20Healthy%20Donor%20-%20No%20Cell%20Sorting%20%2810k%29&lastTouchOfferType=Dataset&product=chromium&redirectUrl=%2Fdatasets%2Fpbmc-from-a-healthy-donor-no-cell-sorting-10-k-1-standard-1-0-0)  

---

## The main analysis includes the following cell types: 

**Cell types retained** :
* B cells
* CD4+ T cells
* CD8+ T cells
* Dendritic cells
* Monocytes
* NK cells
* T cells

**Final dataset after filtering:** A total of 10,989 cells and 578 TF motifs distributed across 7 cell types.
  
---

# cellitac: Model Performance Comparison Across Analytical Frameworks

<p align="center">
    <img src="https://github.com/omicscodeathon/cellitac/blob/main/output/output_Second_Framework_after_dropping/Plots/table02_model_performance.png" alt="cellitac" width="700" />
</p>  
<p align="center">

 
<p align="center">
    <img src="https://github.com/omicscodeathon/cellitac/blob/main/output/output_Second_Framework_after_dropping/Plots/fig07_TF_network_SVM.png" alt="cellitac" width="700" />
</p>  
<p align="center">


---


### Computational Resources

| Pipeline Stage | Hardware Specification |
|----------------|------------------------|
| **Full Pipeline (Preprocessing & ML)** | Personal Laptop: Intel Core 7 240H (16 Cores), 16 GB RAM, 1 TB SSD, NVIDIA GeForce RTX 5060 (8 GB VRAM) (WSL / Linux) |

<br>

### 5. Reproducibility --> need fix

#### Packagies & dependencies :
[all package versions (R - Python) specified for this project](https://github.com/omicscodeathon/cellitac/blob/main/scripts/READme.md)  

----

### 6. License
**License :** [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


## Reporting Issues
To report an issue please use the issues page (https://github.com/omicscodeathon/cellitac/issues). Please check existing issues before submitting a new one.


## Contribute to Project
You can offer to help with the further development of this project by making pull requests on this repo. To do so, fork this repository and make the proposed changes. Once completed and tested, submit a pull request to this repo.

### 7. Contributors  

|   Name   | Affiliation | Role |
|----------|-------------|------|
| **Rana Hamed** | Student, School of Computing and Data Science, Badya University, Cairo, Egypt | Team Lead – Project Management |
| **Syrus Semawule** | African Center of Excellence in Bioinformatics and Data Intensive Sciences, The Infectious Disease Institute, Makerere University, Kampala, Uganda | Bioinformatician – Data Processing & Biological Annotation |
| **Emmanuel Aroma** | Department of Immunology and Molecular Biology, School of Biomedical Sciences, Makerere University, Kampala, Uganda | Bioinformatician – ML Modeling & Pipeline Control |
| **Toheeb Jumah** | Department of Human Anatomy, Faculty of Basic Medical Sciences, College of Medical Sciences, Ahmadu Bello University, Zaria, Nigeria | Bioinformatician – Manuscript Writing & ML Modeling |
| **Olaitan I. Awe** | African Society for Bioinformatics and Computational Biology (ASBCB), Cape Town, South Africa | Project Advisor |

📧 ****Rana Hamed Abu-Zeid**  :** ranahamed2111@gmail.com  
📧 ****Syrus Semawule**  :** semawulesyrus@gmail.com  
📧 ****Emmanuel Aroma**  :** emmatitusaroma@gmail.com  
📧 ****Toheeb Jumah**  :** jumahtoheeb@gmail.com  
📧 ****Olaitan I. Awe, Ph.D.**  :** laitanawe@gmail.com  



---
### Acknowledgments  

We thank the NIH Office of Data Science Strategy for their support before and during the October 2025 Omics Codeathon, co-organized with the African Society for Bioinformatics and Computational Biology (ASBCB).  
We also thank Dr. Awe for his ongoing guidance and all collaborators who contributed to this project.  

---
*This project reflects a collaborative effort towards advancing integrative bioinformatics methods, and we look forward to its continued development and impact within the scientific community.*

