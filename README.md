 # scATAC-tf: 

A reverse TF-centric machine learning framework that classifies peripheral blood mono-nuclear cells (PBMCs) using integrated chromatin accessibility and gene expression data.


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


<p align="center">
  <img src="https://github.com/omicscodeathon/scatactf/blob/main/figures/ChatGPT%20Image%20Sep%2029%2C%202025%2C%2003_39_52%20PM.png" alt="scATAC-tf logo" width="300"  />
</p>

## Table of Contents
1. [Background](#Background)
2. [Presentation Video](#Demonstration_Video)
3. [Code availability](#Code_availability)
4. [License](#License)
5. [Contributors](#Contributors)
<br>

## Background

Single-cell chromatin accessibility sequencing (scATAC-seq) enables genome-wide profiling of regulatory elements at single-cell resolution.Traditional pipelines identify accessible regions first, then infer TF activity, limiting comprehensive understanding of regulatory programs driving cellular identity. This study develops a reverse TF-centric machine learning framework to classify peripheral blood mononuclear cells (PBMCs) using integrated chromatin accessibility and gene expression profiles. Our approach addresses data quality challenges through optimized preprocessing, implements class balancing via SMOTE, and employs ensemble ML methods for robust classification. The resulting computational pipeline enhances single-cell analysis capabilities and provides a systematic approach for discovering TF regulatory networks in immune cell populations.

The current version includes analysis of:

* 5 cell types, namely:

   ★ B-cells
   ★ T-cells
   ★ Monocytes
   ★ NK-cells
   ★ Neutrophils 

This has application in cancer, autoimmune disorders, and developmental conditions since these are not caused by a single broken gene, but by a dysregulation of gene networks.
   
## Presentation Video

<p align="center">
  <a href="https://drive.google.com/file/d/1cNr8JfhEcBRmOS6qTKhtnSOILNbmVBVk/view?usp=sharing">
    <img src="https://github.com/omicscodeathon/scatactf/blob/main/figures/ASBCB-front_image.png" alt="scATAC-tf" width="700" />
  </a>
</p>



### Workflow
<p align="center">
    <img src="https://github.com/omicscodeathon/scatactf/blob/main/figures/ASBCBdrawio_4.png" alt="scATAC-tf" width="700" />
</p>
##  Code Availability  
  

### *Detailed Workflow**  
   👉 [View the Step-by-Step Guide](https://docs.google.com/document/d/1aMClB_2MYsDtn84GWPxdCoQrOLqcOc9K/edit?usp=sharing&ouid=102031592578536141449&rtpof=true&sd=true)  
   A complete walkthrough for running the scATAC-tf pipeline.  

---

###  Demonstration Data  
- Public dataset: [PBMC from a Healthy Donor (10k, 10x Genomics)](https://www.10xgenomics.com/welcome?closeUrl=%2Fdatasets&lastTouchOfferName=PBMC%20from%20a%20Healthy%20Donor%20-%20No%20Cell%20Sorting%20%2810k%29&lastTouchOfferType=Dataset&product=chromium&redirectUrl=%2Fdatasets%2Fpbmc-from-a-healthy-donor-no-cell-sorting-10-k-1-standard-1-0-0)  
- Use this dataset to **test the pipeline** and learn how to correctly organize input data.  

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

## License
**License :** [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



## Reporting Issues
To report an issue please use the issues page (https://github.com/omicscodeathon/scatactf/issues). Please check existing issues before submitting a new one.


## Contribute to Project
You can offer to help with the further development of this project by making pull requests on this repo. To do so, fork this repository and make the proposed changes. Once completed and tested, submit a pull request to this repo.

## Contributors

1. Rana H. Abu-Zeid Student, School of Computing and Data Science, Badya University ,Cairo, Egypt | Team Lead - Project Management   

2. Syrus Semawule, African Center of Excellence in Bioinformatics and Data Intensive Sciences, The Infectious Disease Institute, Makerere University, Kampala, Uganda
 | Bioinformatician - Data processing & Biological Annotation 

3. Emmanuel Aroma, Department of Immunology and Molecular Biology, School of Biomedical Sciences, Makerere University, Kampala, Uganda | Bioinformaticain - ML Modeling & Pipeline Control

4. Toheeb Jumah, Department of Human Anatomy, Faculty of Basic Medical Sciences, College of Medical Sciences, Ahmadu Bello University, Zaria, Nigeria | Bioinformaticain - Manuscript writing & ML Modeling

5. Olaitan I. Awe, African Society for Bioinformatics and Computational Biology (ASBCB), Cape Town, South Africa | Project Advisor
