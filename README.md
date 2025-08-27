<h1 align="center">scatactf</h1>
<h3 align="center">Reverse TF-Centric Modeling of Gene Regulation from scATAC-seq Data</h3>

<h4 align="center"><i>A novel TF-centric framework for analyzing single-cell chromatin accessibility</i></h4>

---

## Overview  

Traditional scATAC-seq analysis follows a forward approach:  
➡️ Open chromatin regions → Identify active TFs → Annotate cell types  

Our **reverse paradigm** starts with known transcription factors (TFs) and predicts their regulatory activity across cell types.  
This enables systematic discovery of TF-driven programs, insights into cellular identity, and mechanisms of disease dysregulation.  

---

## Research Innovation  

**Traditional Workflow:**  
Open Regions → Active TFs → Cell Type Annotation
**Our Reverse Approach:**  
Known TFs → Predict Regulatory Regions → Classify Cell-Specific TF Networks → Discover Novel Mechanisms


---

## Project Scope & Objectives  

- **Primary Goal**: Build machine learning models to classify TFs to their regulatory regions in a cell-type-specific context  
- **Secondary Goal**: Construct comprehensive TF regulatory networks for immune cells  
- **Clinical Impact**: Identify therapeutic targets and biomarkers through TF activity profiling  

---

## Data Sources  
Data Sources & Access
Primary Dataset (Free & Open Access)

10X Genomics Human PBMC Multiome (scATAC-seq + scRNA-seq)

Download link: https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/
Size: ~2-3 GB
Content: 10,847 high-quality immune cells, 156,543 accessible chromatin peaks
Cell types: T-cells, B-cells, Monocytes, NK cells

**Alternative Small Dataset (For Testing)**

3k PBMCs scATAC-seq

Link: https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_pbmc_10k/
Size: ~1 GB
Content: Smaller dataset for initial development

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
Input Features [n_cells × n_features]:
├── TF motif enrichment scores (746 TFs from JASPAR)
├── Chromatin accessibility profiles per cell
├── Gene expression data (from paired scRNA-seq)
├── Cell metadata (cell type, batch, quality metrics)
└── Genomic context features (sequence composition, TSS distance)

Expected Labels:
├── Cell type annotations: ['T_cells', 'B_cells', 'Monocytes', 'NK_cells']
└── TF activity ground truth (from ChIP-seq when available)


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
12. Develop deep learning models (Transformers, GNNs, attention-based)  

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


---

## Computational Framework

- **Programming Languages**: Python (PyTorch, scikit-learn, pandas), R (Seurat, Signac)  
- **Single-cell**: Scanpy, AnnData, episcanpy  $$$$$$$$
- **Genomics**: pybedtools, pyranges, pyfaidx  
- **Visualization**: matplotlib, seaborn, plotly  
- **Infrastructure**: Jupyter Notebook

---

  ## Expected Outputs**

Cell type classification accuracy (>85%)
TF importance rankings per cell type
Basic heatmap visualizations
Model performance metrics

---

## Contributors  

- **[Your Name]** – Principal Investigator | Conceptualization & methodology  
- **[Contributor 2]** – Data Scientist | Model development & validation  
- **[Contributor 3]** – Bioinformatician | Data processing & analysis  
- **[Contributor 4]** – Computational Biologist | Biological interpretation  

Advisors:  
- **[Advisor Name]** – Senior Researcher | Domain expertise  

---

## License  

This project is licensed under the MIT License – see the LICENSE file for details.  

---

## Contact  

- Email: [your-email@institution.edu]  
- GitHub Issues: For bug reports & feature requests  
- Discussions: Join our GitHub discussions tab  

If you find this project useful, please give it a **star** on GitHub! ⭐


