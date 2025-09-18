<p align="center">
  <b>Omics Codeathon General Application - October 2025</b><br>
  Organized by the African Society for Bioinformatics and Computational Biology (ASBCB) with support from the NIH Office of Data Science Strategy.<br>
  Virtual event: October 7–18, 2025 
</p>

<h1 align="center">Reverse TF-Machine Learning Modeling of Gene Regulation from scATAC-seq Data</h1>
<h2 align="center">scATAC-tf</h2>
<h4 align="center"><i>A novel TF-centric framework for analyzing single-cell chromatin accessibility</i></h4>

---
<p align="center">
 <img src="figures/Icon.png" width="600"><br>
</p>
---
<p align="center">
 <img src="figures/ASBCB-front_image.png" width="600"><br>
  <i>Fig.1: Overview of mechanisms influencing chromatin accessibility. Source: <a href="https://www.sc-best-practices.org/chromatin_accessibility/introduction.html">sc-best-practices.org</a> </i>
</p>

---
## Contributors  

- **Rana H. Abu-Zeid** – Team Lead | Computational Biology & Project Management  
- **Syrus Semawule** – Bioinformatician | Data processing & Biological Annotation 
- **Emmanuel Aroma** – Bioinformaticain | ML Modeling & Pipline Control
- **Toheeb Jumah**   – Bioinformaticain | Manuscript writing

Advisors:  
- **Olaitan I. Awe** – Training officer, ASBCB, Cape Town, South Africa. 
---

## Overview

Traditional scATAC-seq analysis follows a forward approach:  

<p align="center">
  <b> Open chromatin regions → Identify active TFs → Annotate cell types</b>
</p>

Our **reverse paradigm** starts with known transcription factors (TFs) and predicts their regulatory activity across cell types.  

<p align="center">
  <b> Known TFs → Predict regulatory regions → Classify cell-specific TF networks → Discover regulatory mechanisms</b>
</p>

This enables systematic discovery of TF-driven programs, insights into cellular identity, and mechanisms of disease dysregulation.

---

## Project Scope & Objectives  

- **Primary Goal**: Build machine learning models to classify TFs to their regulatory regions in a cell-type-specific context  
- **Secondary Goal**: Construct comprehensive TF regulatory networks for immune cells  
- **Clinical Impact**: Identify therapeutic targets and biomarkers through TF activity profiling  

---

## Data Sources  
Data Sources & Access
Primary Dataset (Free & Open Access)
### 📂 **Project Data Overview**

| **Data Type**        | **Source**           | **File Format**                                   |
|----------------------|----------------------|---------------------------------------------------|
| **scATAC-seq**       | 10X Multiome        | `fragments.tsv.gz`, `peaks.bed`                  |
| **scRNA-seq**        | 10X Multiome        | `filtered_feature_bc_matrix.h5`                  |
| **Cell Metadata**    | 10X Multiome        | `singlecell.csv`                                 |
| **TF Motifs**        | JASPAR 2022         | `JASPAR2022_CORE_vertebrates.meme`               |
| **Reference Genome** | UCSC hg38           | `hg38.fa.gz`                                     |
| **Marker Genes**     | PanglaoDB           | `PanglaoDB_markers.tsv`                          |

<p align="center">
  <img src="figures/ASBCB_draw methodology.png" width="600"><br>
  <i>Fig.2: Overview of project methodology</i>
</p>


**Data sets**

10X Genomics Human PBMC Multiome (scATAC-seq + scRNA-seq)

Description: Cryopreserved human peripheral blood mononuclear cells (PBMCs) of a healthy female donor aged 25 were obtained by 10x Genomics from AllCells. Paired ATAC and Gene Expression libraries were generated from the isolated nuclei as described in the Chromium Next GEM Single Cell Multiome ATAC + Gene Expression User Guide (CG000338 Rev A) and sequenced on Illumina Novaseq 6000 v1 Kit (Forward Strand Dual-Index Workflow).

Link to quality checks: [Dataset](https://github.com/omicscodeathon/scatactf/blob/main/data/Study2/pbmc_unsorted_10k%20-%20PBMC%20from%20a%20healthy%20donor%20-%20no%20cell%20sorting%20(10k).html)

- Estimated number of cells: 12,016 
- ATACMean raw read pairs per cell: 38,089 
- ATAC Median high-quality fragments per cell: 9,254 
- ATAC Number of peaks: 63,751 
- GEX Mean raw reads per cell: 69,105 
- GEX Median genes per cell: 1,720 
- GEX Median UMI counts per cell: 3,302 
- Linked genes: 12,576 
- Linked peaks: 46,399 
- Size: ~144 GB

This dataset is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license. 10x citation

**Alternative Small Dataset (For Testing)**
3k PBMCs scATAC-seq
Description: Cryopreserved human peripheral blood mononuclear cells (PBMCs) from a healthy female donor aged 25 were obtained by 10x Genomics from AllCells.

Link to quality check: [Alternative dataset](https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_web_summary.html)

- Estimated number of cells: 3,009
- ATAC Mean raw read pairs per cell: 27,511
- ATAC Median high-quality fragments per cell: 9,877
- ATAC Number of peaks: 81,156
- GEX Mean raw reads per cell: 56,597
- GEX Median genes per cell: 1,494
- GEX Median UMI counts per cell: 2,764
- Linked genes: 6,762
- Linked peaks: 32,613
- Size: ~24 GB

This dataset is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license. 10x citation

**Reference Resources (Open Source)**
JASPAR 2022 – TF motif database
Link: https://jaspar.elixir.no/download/data/2022/CORE/
File: JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt (~50 MB)
Content: 746 TF motifs for vertebrates

Cell Type Markers – From PanglaoDB
Link: https://panglaodb.se/markers.html
Format: CSV file with marker genes per cell type

Human Reference Genome** 
File: hg38.fa.gz (~3 GB)

Link: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
File: hg38.fa.gz (~3 GB)

---
### **Input Features** *(shape: n_cells × n_features)*

- **TF Motif Enrichment Scores:** 746 TFs from JASPAR  
- **Chromatin Accessibility Profiles:** Per-cell accessibility data  
- **Gene Expression Data:** Derived from paired scRNA-seq  
- **Cell Metadata:** Includes cell type, batch information, and quality metrics  
- **Genomic Context Features:** Sequence composition, TSS distance  

---

### **Expected Labels** $$$$$$$

- **Cell Type Annotations:**  
  `['T_cells', 'B_cells', 'Monocytes', 'NK_cells']`  

- **TF Activity Ground Truth:**  
  Based on ChIP-seq data when available  


# Workflow  

### **Block 1: Dataset Collection & Setup**  
1. Download 10X PBMC Multiome dataset (scATAC + scRNA)  
2. Acquire reference databases (JASPAR, ENCODE, GTEx, marker genes)  
3. Configure environment (Python, PyTorch, bedtools, GPU setup)  

---

### **Block 2: Data Preprocessing**  
4. Quality control (remove low-quality cells, doublets, rare peaks)  
5. Normalize accessibility data (log transform, batch correction)  
6. Annotate cell types using RNA marker genes  

---

### **Block 3: Feature Engineering**  
7. Scan DNA sequences for TF motifs (FIMO)  
8. Build TF activity matrix [cells × TFs]  
9. Add genomic features (TSS enrichment, gene expression, sequence context)  

---

### **Block 4: Model Development**  
10. Split dataset into train/validation/test  
11. Train baseline ML models (RF, XGBoost, Logistic Regression)  
12. Develop deep learning models (Transformers, GNNs, attention-based)  $$$$$$ mat not used....

---

### **Block 5: Training & Optimization**  
13. Cross-validation & hyperparameter tuning  
14. Evaluate model performance (Accuracy, F1-score, confusion matrix)  
15. Select best-performing architecture  

---

### **Block 6: TF Activity Prediction**  
16. Generate TF activity scores per cell type  
17. Build TF-cell type regulatory networks  
18. Discover novel TF-region associations  

---

### **Block 7: Biological Validation**  
19. Compare predictions with known TF functions  
20. Validate with external datasets (ENCODE, ChIP-seq)  
21. Perform pathway enrichment & disease relevance analysis  

---

### **Block 8: Biological Interpretation & Annotation**  
22. Identify master regulators via network centrality  
23. Define TF signatures for each immune cell type  
24. Link TF dysregulation to diseases and therapeutic targets  

---

### **Block 9: Results & Visualization**  
25. Generate heatmaps, network graphs, cell type plots  
26. Summarize model performance and biological findings  
27. Package reproducible workflows and deliverables  

---

## Computational Framework

- **Programming Languages**: Python (PyTorch, scikit-learn, pandas), R (Seurat, Signac)  
- **Single-cell**: Scanpy, AnnData, episcanpy  $$$$$$$$
- **Genomics**: pybedtools, pyranges, pyfaidx  
- **Visualization**: matplotlib, seaborn, plotly  
- **Infrastructure**: Jupyter Notebook

---

  ## Expected Outputs**

Cell type classification accuracy 
TF importance rankings per cell type
Basic heatmap visualizations
Model performance metrics

---



## Contact  

- Email: rana.abuzeid@badyau.edu.eg  ,---------,-------------,----------
- GitHub Issues: For bug reports & feature requests  

If you find this project useful, please give it a **STAR** on GitHub! ⭐



