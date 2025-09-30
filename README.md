 # scATAC-tf: 

A reverse TF-centric machine learning framework that classifies peripheral blood mono-nuclear cells (PBMCs) using integrated chromatin accessibility and gene expression data.


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


<p align="center">
  <img src="https://github.com/omicscodeathon/scatactf/blob/main/figures/ChatGPT%20Image%20Sep%2029%2C%202025%2C%2003_39_52%20PM.png" alt="scATAC-tf logo" width="300"  />
</p>

## Table of Contents
1. [Background](#background)
2. [Presentation Video](#2presentation-video)
3. [Workflow](#3workflow)
4. [Code Availability](#4code-avilability)
5. [Reproducibility](#5Reproducibility)
6. [License](#6license)
7. [Contributors](#Contributors)

----
<br>
<p align="center">
  <b>Omics Codeathon General Application - October 2025</b><br>
  Organized by the African Society for Bioinformatics and Computational Biology (ASBCB) with support from the NIH Office of Data Science Strategy.<br>
</p>

----


### 1.Background

Single-cell chromatin accessibility sequencing (scATAC-seq) enables genome-wide profiling of regulatory elements at single-cell resolution.Traditional pipelines identify accessible regions first, then infer TF activity, limiting comprehensive understanding of regulatory programs driving cellular identity. This study develops a reverse TF-centric machine learning framework to classify peripheral blood mononuclear cells (PBMCs) using integrated chromatin accessibility and gene expression profiles. Our approach addresses data quality challenges through optimized preprocessing, implements class balancing via SMOTE, and employs ensemble ML methods for robust classification. The resulting computational pipeline enhances single-cell analysis capabilities and provides a systematic approach for discovering TF regulatory networks in immune cell populations.

----

 
### 2.Presentation Video

<p align="center">
  <a href="https://drive.google.com/file/d/1cNr8JfhEcBRmOS6qTKhtnSOILNbmVBVk/view?usp=sharing">
    <img src="https://github.com/omicscodeathon/scatactf/blob/main/figures/ASBCB-front_image.png" alt="scATAC-tf" width="700" />
  </a>
</p>



### 3.Workflow  

<p align="center">
    <img src="https://github.com/omicscodeathon/scatactf/blob/main/figures/ASBCBdrawio_4.png" alt="scATAC-tf" width="700" />
</p>  
<p align="center">
    <b>Figure 1.</b> Workflow of the methods employed in this study
</p>

 
  

### *Detailed Workflow*
   👉 [View the Step-by-Step Guide](https://docs.google.com/document/d/1aMClB_2MYsDtn84GWPxdCoQrOLqcOc9K/edit?usp=sharing&ouid=102031592578536141449&rtpof=true&sd=true)  
   A complete walkthrough for running the scATAC-tf pipeline.  

## Pipeline Architecture

```mermaid
graph LR
    A[Raw 10X Data] --> B[Quality Control]
    B --> C[Feature Engineering & integration]
    C --> D[ML Building]
    D --> E[Evalutaion Metrics]
    E --> F[Validation & interpretation]
---

### 4.Code Avilability:

All scripts for the **scATAC-tf** project (Python & R) are available in the repository:

👉 Browse the scripts: [scripts/](scripts/)

---

###  Demonstration Data  
- Public dataset: [PBMC from a Healthy Donor (10k, 10x Genomics)](https://www.10xgenomics.com/welcome?closeUrl=%2Fdatasets&lastTouchOfferName=PBMC%20from%20a%20Healthy%20Donor%20-%20No%20Cell%20Sorting%20%2810k%29&lastTouchOfferType=Dataset&product=chromium&redirectUrl=%2Fdatasets%2Fpbmc-from-a-healthy-donor-no-cell-sorting-10-k-1-standard-1-0-0)  

---

## The main analysis includes the following cell types: 

**Cell types retained** :
* B cells 
* Monocytes 
* NK cells
* T cells 

**Excluded rare cell types** (<10 samples):
* HSC-G-CSF 
* Pre-B cells CD34- 

**Final dataset after filtering:** about ≈ 1,400 cells across 4 cell types
  
---

###  Computational Framework  

| Language | Key Packages |
|----------|--------------|
| **Python** | PyTorch, scikit-learn, pandas |
| **R**      | Seurat, Signac |

---

###  Computational Resources  

| Step | Recommended Resources |
|------|------------------------|
| **Pre-processing** | 3.6 GHz 10-Core Intel Core i9, 64 GB RAM, 10 GB storage |
| **Modeling & Scripts**     | Open-source code available (optimized for scalability) |


----

<br>

### 5.Reproducibility

#### Random Seeds
All scripts use fixed random seed (42) for reproducibility
#### Packagies & dependencies
[all package versions (R - Python) specified for this project](https://docs.google.com/document/d/1aMClB_2MYsDtn84GWPxdCoQrOLqcOc9K/edit?usp=sharing&ouid=102031592578536141449&rtpof=true&sd=true)  


----

Results of core Framework_1 :
<p align="center">

| Model          | Accuracy   | Precision  | Recall     | F1 Score   | AUC        |
|----------------|------------|------------|------------|------------|------------|
| Random Forest  | 0.9725     | 0.9718     | 0.9725     | 0.9721     | 0.9983     |
| XGBoost        | 0.9828     | 0.9827     | 0.9828     | 0.9825     | 0.9992     |
| SVM            | 0.9794     | 0.9800     | 0.9794     | 0.9780     | 0.9977     |

</p>


### 6.License
**License :** [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



## Reporting Issues
To report an issue please use the issues page (https://github.com/omicscodeathon/scatactf/issues). Please check existing issues before submitting a new one.


## Contribute to Project
You can offer to help with the further development of this project by making pull requests on this repo. To do so, fork this repository and make the proposed changes. Once completed and tested, submit a pull request to this repo.

##  Contributors  

|   Name   | Affiliation | Role |
|----------|-------------|------|
| **Rana Hamed** | Student, School of Computing and Data Science, Badya University, Cairo, Egypt | Team Lead – Project Management |
| **Syrus Semawule** | African Center of Excellence in Bioinformatics and Data Intensive Sciences, The Infectious Disease Institute, Makerere University, Kampala, Uganda | Bioinformatician – Data Processing & Biological Annotation |
| **Emmanuel Aroma** | Department of Immunology and Molecular Biology, School of Biomedical Sciences, Makerere University, Kampala, Uganda | Bioinformatician – ML Modeling & Pipeline Control |
| **Toheeb Jumah** | Department of Human Anatomy, Faculty of Basic Medical Sciences, College of Medical Sciences, Ahmadu Bello University, Zaria, Nigeria | Bioinformatician – Manuscript Writing & ML Modeling |
| **Olaitan I. Awe** | African Society for Bioinformatics and Computational Biology (ASBCB), Cape Town, South Africa | Project Advisor |

📧 ****Rana Hamed Abu-Zeid**  :** rana.hamed@badyauni.edu.eg  
📧 ****Syrus Semawule**  :** semawulesyrus@gmail.com  
📧 ****Emmanuel Aroma**  :** emmatitusaroma@gmail.com  
📧 ****Toheeb Jumah**  :** jumahtoheeb@gmail.com  


---
### Acknowledgments  

We thank the NIH Office of Data Science Strategy for their support before and during the April 2025 Omics Codeathon, co-organized with the African Society for Bioinformatics and Computational Biology (ASBCB).  
We also thank Dr. Awe for his ongoing guidance and all collaborators who contributed to this project.  

---
*This project reflects a collaborative effort towards advancing integrative bioinformatics methods, and we look forward to its continued development and impact within the scientific community.*

