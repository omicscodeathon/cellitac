
## Project Overview

**Core Objective**: Build machine learning models to predict transcription factor (TF) regulatory activity from single-cell chromatin accessibility data (scATAC-seq) in immune cells.

**Innovation**: Reverse paradigm approach - Start with known TFs → Predict their regulatory activity → Classify cell-specific TF networks

**Expected Impact**: Identify therapeutic targets and biomarkers through TF activity profiling in immune cells.

**Timeline**: 3-4 weeks (15-18 working days with breaks every 3 days)

---

## Module 1: Data Acquisition & Setup

### **Input Files (Download These)**:
```
📁 Raw Data Sources/
├── 10X_PBMC_Multiome/
│   ├── fragments.tsv.gz           # scATAC fragments (3-4 GB)
│   ├── filtered_feature_bc_matrix.h5  # scRNA + ATAC matrix (1-2 GB)
│   ├── singlecell.csv             # Cell metadata (10 MB)
│   └── peaks.bed                  # Called peaks (5 MB)
├── JASPAR2022_CORE_vertebrates.meme    # TF motifs (50 MB)
├── hg38.fa.gz                     # Reference genome (3 GB)
└── PanglaoDB_markers.tsv          # Cell type markers (5 MB)
```

**Download Links**:
- 10X Data: https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/
- JASPAR: https://jaspar.elixir.no/download/data/2022/CORE/
- hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
- PanglaoDB: https://panglaodb.se/markers.html

### **What We Do**:
1. Download and organize all required datasets
2. Set up R environment (ArchR, Seurat, Signac, chromVAR)
3. Set up Python environment (scanpy, muon, scikit-learn, XGBoost)
4. Verify data integrity and file formats

### **Tools**: R + Python setup

### **Output Files**:
```
📁 01_Raw_Data/
├── fragments.tsv.gz
├── filtered_feature_bc_matrix.h5
├── singlecell.csv
├── peaks.bed
├── jaspar_motifs.meme
├── genome_hg38.fa
├── marker_genes.tsv
└── data_inventory.txt          # File verification log
```

**Time Required**: 2-3 days

---

## Module 2: Data Preprocessing & Quality Control

### **Input Files**:
```
📁 01_Raw_Data/
├── fragments.tsv.gz        # MAIN INPUT: Process THIS for accessibility
├── filtered_feature_bc_matrix.h5  # MAIN INPUT: Process THIS for RNA + cell annotation
├── singlecell.csv          # Use THIS for initial cell metadata
└── marker_genes.tsv        # Use THIS for cell type annotation
```

### **What We Do**:

#### **Step 2.1: Cell Quality Control**
**Process**: `fragments.tsv.gz` + `singlecell.csv`
- Filter cells: >2000 fragments, TSS enrichment >6
- Remove doublets using ArchR's doublet detection
- Calculate Fraction of Reads in Peaks (FRiP) >0.3

#### **Step 2.2: Peak Quality Control**  
**Process**: `peaks.bed` + `fragments.tsv.gz`
- Remove ENCODE blacklist regions
- Filter peaks present in <1% of cells
- Call peaks using MACS2 (pseudo-bulk approach)

#### **Step 2.3: Cell Type Annotation**
**Process**: `filtered_feature_bc_matrix.h5` + `marker_genes.tsv`
- Extract scRNA-seq data from matrix.h5
- Apply reference-based annotation using PanglaoDB markers
- Target cell types: T-cells, B-cells, Monocytes, NK cells

### **Primary Tool**: R (ArchR/Seurat)

### **R Code Example**:
```r
# Load main accessibility data
arrows <- createArrowFiles("fragments.tsv.gz")
proj <- ArchRProject(ArrowFiles = arrows)

# Quality control
proj <- filterDoublets(proj)
proj <- addIterativeLSI(proj)

# Cell type annotation using RNA
proj <- addGeneIntegrationMatrix(proj, useMatrix = "GeneScoreMatrix")
```

### **Output Files**:
```
📁 02_Processed_Data/
├── clean_atac_matrix.h5      # Filtered accessibility matrix [cells × peaks]
├── clean_rna_matrix.h5       # Filtered RNA matrix [cells × genes]
├── cell_metadata_annotated.csv  # QC metrics + cell type labels
├── filtered_peaks.bed           # High-quality peaks only
├── quality_control_report.html  # QC statistics and plots
└── cell_type_distribution.csv   # Cell count per type
```

**Time Required**: 3-4 days

---

## Module 3: Feature Engineering

### **Input Files**:
```
📁 02_Processed_Data/
├── clean_atac_matrix.h5      # MAIN INPUT: Normalize THIS
├── filtered_peaks.bed        # MAIN INPUT: Scan THIS for motifs
├── cell_metadata_annotated.csv  # Use THIS for cell stratification

📁 01_Raw_Data/
├── jaspar_motifs.meme        # MAIN INPUT: Use THIS for motif scanning
└── genome_hg38.fa           # Use THIS for sequence extraction
```

### **What We Do**:

#### **Step 3.1: Chromatin Accessibility Features**
**Process**: `clean_atac_matrix.h5`
- Apply TF-IDF normalization to accessibility counts
- Perform LSI dimensionality reduction (50 components)
- Remove depth-correlated first component

#### **Step 3.2: TF Motif Features** 
**Process**: `filtered_peaks.bed` + `jaspar_motifs.meme` + `genome_hg38.fa`
- Extract DNA sequences for each peak from reference genome
- Scan sequences for JASPAR motifs using FIMO
- Calculate motif enrichment scores per cell (chromVAR approach)
- Generate 746-dimensional TF motif activity matrix

#### **Step 3.3: Genomic Context Features**
**Process**: `filtered_peaks.bed` + `genome_hg38.fa`
- Calculate distance to nearest TSS
- Compute peak width and signal strength
- Extract sequence composition (GC content, k-mer frequencies)

#### **Step 3.4: Integration Features**
**Process**: `clean_atac_matrix.h5` + `clean_rna_matrix.h5`
- Generate gene activity scores (peak-to-gene linking)
- Include RNA expression of TF genes from paired data

### **Primary Tool**: R (ArchR + chromVAR) → Export to Python

### **R Code Example**:
```r
# Main motif processing
motifPositions <- getPositions(proj, name = "Motif")
proj <- addMotifAnnotations(proj, motifSet = "jaspar2022")
proj <- addDeviationsMatrix(proj, peakAnnotation = "Motif")
```

### **Output Files**:
```
📁 03_Features/
├── tfidf_accessibility_matrix.csv    # Normalized accessibility [cells × peaks]
├── tf_motif_scores_746.csv          # TF motif scores [cells × 746 TFs]
├── genomic_context_features.csv      # Context features [cells × genomic_features]
├── gene_activity_scores.csv         # Peak-to-gene scores [cells × genes]
├── combined_feature_matrix.csv       # ALL features combined [cells × all_features]
├── feature_descriptions.txt          # Description of each feature column
└── motif_annotation_summary.csv     # Motif scanning statistics
```

**Time Required**: 4-5 days

---

## Module 4: Machine Learning Model Development

### **Input Files**:
```
📁 03_Features/
├── combined_feature_matrix.csv   # MAIN INPUT: Train models on THIS
├── tf_motif_scores_746.csv      # Alternative input: TF-focused features
└── cell_metadata_annotated.csv   # MAIN INPUT: Use THIS for target labels
```

### **What We Do**:

#### **Step 4.1: Data Preparation**
**Process**: `combined_feature_matrix.csv` + `cell_metadata_annotated.csv`
- Split data into train/validation/test (70/15/15)
- Stratify by cell type to ensure balanced representation
- Handle missing values and feature scaling

#### **Step 4.2: Baseline Models**
**Process**: Training features → Simple models
- Logistic Regression with L1/L2 regularization
- Random Forest classifier
- Purpose: Establish baseline performance and feature importance

#### **Step 4.3: Advanced Models**
**Process**: Same training features → Complex models
- XGBoost/LightGBM (handle non-linear interactions)
- Support Vector Machines with RBF kernel
- Purpose: Capture complex TF-chromatin relationships

#### **Step 4.4: Model Optimization**
**Process**: All models → Hyperparameter tuning
- 5-fold cross-validation for model selection
- Grid search for optimal parameters
- Early stopping to prevent overfitting

### **Primary Tool**: Python (scikit-learn, XGBoost)

### **Python Code Example**:
```python
# Main processing
features = pd.read_csv("combined_feature_matrix.csv")
labels = pd.read_csv("cell_metadata_annotated.csv")['cell_type']

# Train models
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier

rf_model = RandomForestClassifier().fit(X_train, y_train)
xgb_model = XGBClassifier().fit(X_train, y_train)
```

### **Target Variables**:
- Cell type classification (4 classes: T-cells, B-cells, Monocytes, NK cells)
- TF activity scores (continuous regression per TF)
- Active vs. inactive TF states (binary classification per TF)

### **Output Files**:
```
📁 04_Models/
├── logistic_regression_model.pkl    # Trained baseline model
├── random_forest_model.pkl          # Trained RF model  
├── xgboost_model.pkl               # Trained XGB model (likely best)
├── svm_model.pkl                   # Trained SVM model
├── model_parameters.json            # All hyperparameters used
├── cross_validation_results.csv     # CV performance per model
├── feature_importance_all_models.csv # Feature rankings per model
└── training_logs.txt               # Training progress and timing
```

**Time Required**: 3-4 days

---

## Module 5: Model Evaluation & Validation

### **Input Files**:
```
📁 04_Models/
├── logistic_regression_model.pkl   # Load THESE for predictions
├── random_forest_model.pkl
├── xgboost_model.pkl
└── svm_model.pkl

📁 03_Features/
├── combined_feature_matrix.csv     # Use THIS for test predictions
└── cell_metadata_annotated.csv     # Use THIS for true labels (ground truth)

📁 01_Raw_Data/
└── marker_genes.tsv               # Use THIS for biological validation
```

### **What We Do**:

#### **Step 5.1: Statistical Performance**
**Process**: All trained models + test data
- Generate predictions on held-out test set
- Calculate classification metrics: Accuracy, F1-score, Precision, Recall
- Compute AUROC and AUPRC (robust to class imbalance)
- Create confusion matrices per cell type

#### **Step 5.2: Clustering Quality Assessment**
**Process**: Model predictions vs. true cell types
- Adjusted Rand Index (ARI) - primary metric for cluster comparison
- Adjusted Mutual Information (AMI)
- Silhouette scores for cluster separation

#### **Step 5.3: Biological Validation**
**Process**: TF predictions + known biology from `marker_genes.tsv`
- Test recovery of established TF-cell type associations
- Examples: SPI1 in monocytes, PAX5 in B-cells, CD3 in T-cells
- Compare with external ChIP-seq data (if available)

#### **Step 5.4: Model Comparison**
**Process**: All model outputs → Select best performer
- Rank models by ARI and biological validation scores
- Select final model for downstream interpretation

### **Primary Tool**: Python (scikit-learn metrics) + R (visualization)

### **Python Code Example**:
```python
# Main evaluation processing
test_features = pd.read_csv("combined_feature_matrix.csv").iloc[test_indices]
true_labels = pd.read_csv("cell_metadata_annotated.csv").iloc[test_indices]['cell_type']

# Load and evaluate models
import pickle
from sklearn.metrics import accuracy_score, adjusted_rand_score

xgb_model = pickle.load(open("xgboost_model.pkl", "rb"))
predictions = xgb_model.predict(test_features)
accuracy = accuracy_score(true_labels, predictions)
ari = adjusted_rand_score(true_labels, predictions)
```

### **Performance Targets**:
- Cell type classification: >85% accuracy
- TF activity prediction: ARI >0.7 with known patterns
- Novel TF associations: >70% validation rate

### **Output Files**:
```
📁 05_Evaluation/
├── model_predictions_all.csv         # Predictions from all models
├── performance_metrics_summary.csv   # Accuracy, F1, AUROC, ARI for each model
├── confusion_matrices.png           # Visual confusion matrices
├── roc_curves_all_models.png        # ROC curves comparison
├── feature_importance_ranked.csv     # Top features driving predictions
├── biological_validation_results.csv # Known TF associations recovered
├── model_comparison_report.html      # Complete evaluation summary
└── best_model_selection.txt         # Final model choice with justification
```

**Time Required**: 2-3 days

---

## Module 6: Biological Interpretation & Visualization

### **Input Files**:
```
📁 05_Evaluation/
├── model_predictions_all.csv       # Use THESE for TF activity analysis
├── feature_importance_ranked.csv   # Use THIS to identify key TFs
├── biological_validation_results.csv # Use THIS for known patterns
└── best_model_selection.txt        # Use THIS to focus on best model

📁 03_Features/
├── tf_motif_scores_746.csv     # MAIN INPUT: Use THIS for TF network construction
├── cell_metadata_annotated.csv # Use THIS for cell type stratification
└── combined_feature_matrix.csv  # Use THIS for SHAP analysis

📁 04_Models/
└── xgboost_model.pkl          # Use THIS (best model) for interpretation
```

### **What We Do**:

#### **Step 6.1: TF Network Construction**
**Process**: `tf_motif_scores_746.csv` + `cell_metadata_annotated.csv`
- Calculate cell-type-specific TF activity profiles
- Identify master regulator TFs per cell type
- Build TF-TF interaction networks
- Map TF regulatory hierarchies

#### **Step 6.2: Feature Importance Analysis**
**Process**: `best_model.pkl` + `combined_feature_matrix.csv`
- Compute SHAP values for model interpretability
- Identify key chromatin regions driving TF activity
- Rank motifs by predictive importance across cell types
- Perform in-silico perturbation analysis

#### **Step 6.3: Pathway Integration**
**Process**: TF activity profiles + pathway databases
- Connect TF activity to biological pathways
- Disease association analysis using literature
- Therapeutic target identification and prioritization

#### **Step 6.4: Visualization Generation**
**Process**: All results → Create publication-ready figures
- TF activity heatmaps per cell type
- Regulatory network graphs
- UMAP/tSNE plots colored by TF activity
- Feature importance plots and summary statistics

### **Tools**: R (ggplot2, ComplexHeatmap) + Python (SHAP, matplotlib)

### **R + Python Code Example**:
```r
# R: TF network analysis
tf_activity <- read.csv("tf_motif_scores_746.csv")
cell_types <- read.csv("cell_metadata_annotated.csv")

# Create cell-type specific TF profiles
tf_profiles <- aggregate(tf_activity, by=list(cell_types$cell_type), mean)
```

```python
# Python: SHAP interpretation
import shap
import pickle

model = pickle.load(open("xgboost_model.pkl", "rb"))
explainer = shap.Explainer(model)
shap_values = explainer(test_features)
```

### **Output Files**:
```
📁 06_Results/
├── tf_activity_heatmap.png          # MAIN RESULT: TF activity per cell type
├── regulatory_networks.png          # TF interaction network graphs
├── cell_type_tf_profiles.csv        # TF activity scores per cell type
├── master_regulators_ranked.csv     # Top TFs per cell type with scores
├── shap_feature_importance.png      # Model interpretability plots
├── novel_tf_associations.csv        # New TF-cell discoveries
├── therapeutic_targets_prioritized.csv # Clinical relevance ranking
├── pathway_enrichment_results.csv   # TF-pathway associations
├── umap_tf_activity.png            # Dimensionality reduction plots
└── final_comprehensive_report.html  # Complete analysis summary
```

**Time Required**: 3-4 days

---

## Complete File Flow Summary

### **Critical Data Transformation Path**:
```
fragments.tsv.gz (3-4 GB) 
    ↓ Module 2: QC & filtering
clean_atac_matrix.h5 (~1 GB)
    ↓ Module 3: TF-IDF + motif scanning  
tf_motif_scores_746.csv (~50 MB)
    ↓ Module 4: ML training
trained_models.pkl (~10 MB)
    ↓ Module 5: Predictions
model_predictions.csv (~5 MB)
    ↓ Module 6: Interpretation
tf_activity_profiles.csv (~1 MB) → FINAL BIOLOGICAL INSIGHTS
```

### **Most Critical Files to Monitor**:
1. **fragments.tsv.gz** → Raw accessibility data (largest file)
2. **filtered_feature_bc_matrix.h5** → RNA data for cell annotation
3. **tf_motif_scores_746.csv** → Core features for ML models
4. **combined_feature_matrix.csv** → Final training input
5. **tf_activity_heatmap.png** → Main result visualization

---

## Implementation Strategy

### **Phase 1: R-based Processing (Week 1)**
**Days 1-3**: Modules 1-2 (Setup + Preprocessing)
- Focus on `fragments.tsv.gz` and `matrix.h5` processing
- Generate clean, annotated datasets

**Day 4**: Break

**Days 5-7**: Module 3 (Feature Engineering)  
- Focus on motif scanning and feature matrix creation

### **Phase 2: Python-based Modeling (Week 2)**
**Days 1-3**: Module 4 (ML Development)
- Train models on `combined_feature_matrix.csv`
- Optimize hyperparameters

**Day 4**: Break

**Days 5-7**: Module 5 (Evaluation)
- Comprehensive model evaluation and validation

### **Phase 3: Integration & Results (Week 3)**
**Days 1-3**: Module 6 (Interpretation)
- Biological analysis and visualization

**Day 4**: Break

**Days 5-7**: Final documentation and results packaging

---

## Success Metrics

**Technical Success**:
- All modules complete without errors
- Models achieve >85% cell type classification accuracy
- Feature matrices properly constructed from raw files

**Biological Success**:
- Recover known TF-cell type associations from literature
- Generate interpretable TF regulatory networks
- Identify actionable therapeutic targets

**File Management Success**:
- Clear data provenance from raw to final results
- Reproducible workflow with documented file transformations
- Organized output structure for easy sharing and publication

This unified approach provides a clear roadmap with explicit file dependencies, making it straightforward to implement and debug each step while maintaining focus on your core research objectives.
