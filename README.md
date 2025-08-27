**Reverse TF-Centric Modeling of Gene Regulation from scATAC-seq**

---

##  Overview  
Traditional single-cell chromatin accessibility analysis follows a forward approach:  
➡identify open chromatin regions → find transcription factors → annotate cell types.  

Our research introduces a **reverse classification paradigm** that:  
 starts with known transcription factors (TFs) and predicts their cell-type-specific regulatory activity patterns.  

This innovative approach enables:  
- systematic discovery of **TF-driven regulatory programs**  
- new insights into **cellular identity determination** and **disease mechanisms**  

---
###  Traditional Workflow  
Open Regions → Active TFs → Cell Type Annotation

##  Research Innovation  :
Known TFs → Regulatory Region Prediction → Cell-Specific TF Networks → Novel Regulatory Discovery

---

##  Project Scope & Objectives  

- **Primary Goal**: Develop machine learning models to classify transcription factors to their corresponding regulatory regions in cell-type-specific contexts  
- **Secondary Goal**: Create comprehensive TF regulatory network maps for human immune cells  
- **Clinical Impact**: Identify therapeutic targets and biomarkers through TF activity profiling  

---

## Data Sources  

### 🔹 Primary Datasets  
- **10X Genomics Human PBMC Multiome** (scATAC-seq + scRNA-seq)  
  - 10,847 high-quality immune cells  
  - 156,543 accessible chromatin peaks  
  - 4 major cell types: *T-cells, B-cells, Monocytes, NK cells*  

### 🔹 Reference Resources  
- **JASPAR 2022** – Transcription factor binding motif database  
- **ENCODE Project** – Chromatin accessibility & histone modification data  
- **GTEx Atlas** – Tissue-specific gene expression profiles  
- **ChIP-Atlas** – TF binding site validation data  

---

## Input & Output Architecture  

### 🔹 Input Features  
**Multi-modal Feature Matrix** `[n_cells × n_features]`  
