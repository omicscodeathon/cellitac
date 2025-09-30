# Workflow  

<p align="center">
    <img src="https://github.com/omicscodeathon/scatactf/blob/main/figures/ASBCBdrawio_4.png" alt="scATAC-tf" width="700" />
</p>
--- 


### **Block 1: Dataset Collection & Setup**  
1. Download 10X PBMC Multiome dataset  
2. Configure environment (Python, RScript)  

---

### **Block 2: Preprocessing & QC**  
1.Cell Filtering Criteria
2.Gene Filtering
3.Quality control metrics
4.Peak Filtering
5.TF-IDF normalization (ATAC)
6.Automated Cell Type Annotation

---

### **Block 3: Feature Engineering**  
***1.RNA Features:***
    1.1.Select top variable genes
    1.2.Log-normalization
    
***2.ATAC Features:***
    3.1.Select top accessible peaks
    3.2.Normalization with LSI
    3.3.Multi-modal integration
    3.4.Feature selection & scaling

---

### **Block 4: Model Development**  
1. Split dataset into train/validation/test
2. Train baseline ML models
   
---

### **Block 5: Evaluation**  
1. Cross-validation & hyperparameter tuning  
1. Evaluate model performance (Accuracy, F1-score, confusion matrix)  

---

### **Block 6: validation & interpretation**  
1.Feature importance analysis
2.TF-cell type networks
3.Visualization & documentation

---
