# Complete scATAC-tf Analysis Pipeline
## From Data Acquisition to Advanced Visualization

---

## Overview

This document details a comprehensive computational pipeline for multi-modal single-cell RNA-seq and ATAC-seq data integration, cell-type classification, and advanced biological interpretation.

The workflow combines **R-based preprocessing** with **Python-based machine learning** and visualization.

---

## Table of Contents

- [Phase 1: Data Acquisition and Preprocessing (R Scripts)](#phase-1-data-acquisition-and-preprocessing-r-scripts)
  - [Step 1: Data Download and Setup](#step-1-data-download-and-setup)
  - [Step 2: RNA Data Processing](#step-2-rna-data-processing-r-file-part-1)
  - [Step 3: ATAC Data Processing](#step-3-atac-data-processing-r-file-part-2)
  - [Step 4: Multi-Modal Data Integration](#step-4-multi-modal-data-integration)
- [Phase 2: Machine Learning and Classification (Python Script)](#phase-2-machine-learning-and-classification-python-script)
  - [Step 5: Environment Setup and Data Loading](#step-5-environment-setup-and-data-loading)
  - [Step 6: Class Imbalance Analysis and Handling](#step-6-class-imbalance-analysis-and-handling)
  - [Step 7: Feature Engineering and Selection](#step-7-feature-engineering-and-selection)
  - [Step 8: Model Training and Evaluation](#step-8-model-training-and-evaluation)
  - [Step 9: Cross-Validation and Overfitting Analysis](#step-9-cross-validation-and-overfitting-analysis)
- [Phase 3: Advanced Analysis and Visualization](#phase-3-advanced-analysis-and-visualization)
  - [Step 10: Feature Importance Analysis](#step-10-feature-importance-analysis)
  - [Step 11: Cell Type-Specific Feature Analysis](#step-11-cell-type-specific-feature-analysis)
  - [Step 12: Visualization Suite Creation](#step-12-visualization-suite-creation)
  - [Step 13: Results Export and Documentation](#step-13-results-export-and-documentation)
- [Summary of Critical Parameters](#summary-of-critical-parameters-for-publication)
- [Contact](#contact)

---

## PHASE 1: Data Acquisition and Preprocessing (R Scripts)

### Step 1: Data Download and Setup

**Script:** Initial setup and data acquisition  
**Packages Used:** Base R functions

#### Parameters:
- **Data source:** 10x Genomics PBMC dataset (10k cells)
- **Required files:** 
  - `filtered_feature_bc_matrix.h5`
  - `per_barcode_summary_metrics.csv`
  - `atac_fragments.tsv.gz`
  - `atac_peaks.bed`

#### Results:
- Multi-modal dataset with both RNA and ATAC components available
- Raw data loaded and ready for quality control

---

### Step 2: RNA Data Processing (R file Part 1)

**Script:** Team 1 RNA processing pipeline  
**Packages Used:**
- Seurat v5.3.0
- hdf5r
- dplyr
- ggplot2
- SingleR v2.0.0
- celldex

#### Key Parameters:
- **Minimum features per cell:** 500 genes
- **Maximum features per cell:** 7000 genes (doublet removal)
- **Minimum UMI counts:** 1000
- **Maximum mitochondrial percentage:** 20%
- **Variable features:** 2000 most variable genes
- **PCA dimensions:** 1-15
- **Clustering resolution:** 0.5
- **Normalization:** LogNormalize with scale factor 10,000

#### Processing Steps:

**1. Quality Control Filtering:**
- Remove cells with <500 or >7000 genes
- Remove cells with >20% mitochondrial genes
- Remove genes expressed in <5 cells

**2. Normalization and Scaling:**
- Log-normalization with scale factor 10,000
- Variance-stabilizing transformation for feature selection
- Z-score scaling for PCA

**3. Dimensionality Reduction:**
- PCA on top 2000 variable genes
- UMAP visualization using first 15 PCs
- t-SNE for additional visualization

**4. Cell Type Annotation:**
- Automated annotation using SingleR
- Reference: HumanPrimaryCellAtlasData
- Manual validation of cluster markers

#### Results:
- Processed RNA data: X cells × 2000 variable genes
- Cell type annotations for all cells
- Quality control metrics and visualizations
- Dimensionality reduction coordinates

---

### Step 3: ATAC Data Processing (R file part 2)

**Script:** Team 2 ATAC processing pipeline  
**Packages Used:**
- Signac v1.15.0
- GenomicRanges
- EnsDb.Hsapiens.v75
- biovizBase

#### Key Parameters:
- **Minimum peak counts per cell:** 1000
- **Maximum peak counts per cell:** 100,000
- **Minimum TSS enrichment:** 2
- **Maximum nucleosome signal:** 4
- **Minimum reads in peaks:** 15%
- **Maximum blacklist ratio:** 0.05
- **Top accessible peaks:** 5000

#### Processing Steps:

**1. Quality Control Filtering:**
- Filter cells based on fragment counts and quality metrics
- Remove low-quality peaks and blacklist regions
- Calculate TSS enrichment and nucleosome banding scores

**2. Normalization:**
- TF-IDF normalization for accessibility data
- Feature selection of top variable peaks

**3. Dimensionality Reduction:**
- Latent Semantic Indexing (LSI) for ATAC data
- UMAP visualization using LSI dimensions 2-30

#### Results:
- Processed ATAC data: X cells × 5000 accessible peaks
- Quality control metrics
- Chromatin accessibility profiles
- Peak accessibility scores

---

### Step 4: Multi-Modal Data Integration

**Script:** ASBCS_step2_Syrus.r  
**Packages Used:** Base R, data manipulation functions

#### Key Parameters:
- **Train/test split ratio:** 70/30
- **Stratified sampling by cell type**
- **Random seed:** 42

#### Processing Steps:

**1. Cell Alignment:**
- Identify common cells between RNA and ATAC modalities
- Create unified feature matrix combining both modalities

**2. Train/Test Splitting:**
- Stratified sampling maintaining cell type proportions
- Separate training and testing sets

**3. Data Export:**
- Combined feature matrix (cells × features)
- Cell type labels and metadata
- Train/test split assignments
- Feature type annotations (RNA vs ATAC)

#### Results:
- Combined multi-modal dataset: X cells × 7000 features (2000 RNA + 5000 ATAC)
- Stratified train/test splits with balanced cell type representation
- Python-ready CSV files for downstream machine learning

---

## PHASE 2: Machine Learning and Classification (Python Script)

### Step 5: Environment Setup and Data Loading

**Script:** Preliminary_result_v2.py  
**Packages Used:**
- pandas v2.2.3
- numpy v2.1.2
- scikit-learn v1.6.1
- matplotlib v3.10.1
- seaborn v0.13.2
- imblearn (imbalanced-learn v0.14.0)

#### Processing Steps:
1. Load multi-modal feature matrix
2. Load cell type labels and train/test splits
3. Data validation and integrity checking

#### Results:
- Validated dataset ready for machine learning
- Class distribution analysis

---

### Step 6: Class Imbalance Analysis and Handling

**Packages Used:** imblearn, matplotlib, seaborn

#### Key Parameters:
- **Imbalance threshold:** 2.0 (max/min class ratio)
- **SMOTE variant:** BorderlineSMOTE
- **SMOTE parameters:** k_neighbors=3, kind="borderline-1"

#### Processing Steps:

**1. Imbalance Detection:**
- Calculate class distribution and imbalance ratios
- Determine if balancing is required

**2. SMOTE Application (if needed):**
- Generate synthetic samples for minority classes
- Maintain feature space integrity

#### Results:
- Balanced training dataset (if balancing applied)
- Class distribution visualizations
- Imbalance analysis metrics

---

### Step 7: Feature Engineering and Selection

**Packages Used:** scikit-learn preprocessing and feature selection

#### Key Parameters:
- **Correlation threshold:** 0.95 (for redundant feature removal)
- **Feature selection:** SelectKBest with f_classif
- **Selected features:** min(1000, total_features/2)
- **Feature scaling:** StandardScaler

#### Processing Steps:

**1. Feature Cleaning:**
- Remove zero-variance features
- Remove highly correlated features (>0.95)

**2. Statistical Feature Selection:**
- Univariate feature selection using F-statistics
- Select top discriminative features

**3. Feature Scaling:**
- Z-score normalization for algorithms requiring scaled input

#### Results:
- Reduced feature set optimized for classification
- Scaled features for algorithm compatibility

---

### Step 8: Model Training and Evaluation

**Packages Used:** scikit-learn ensemble, linear_model, svm, neural_network, xgboost

#### Key Parameters:

**Random Forest:**
- n_estimators: 100
- max_depth: 10
- min_samples_split: 20
- min_samples_leaf: 10
- class_weight: 'balanced'

**XGBoost:**
- n_estimators: 100
- max_depth: 6
- learning_rate: 0.1
- subsample: 0.8
- reg_alpha: 0.1, reg_lambda: 1.0

**SVM:**
- kernel: 'rbf'
- C: 1.0
- class_weight: 'balanced'

#### Processing Steps:

**1. Model Training:**
- Train multiple algorithms with regularization
- Apply appropriate data scaling per algorithm

**2. Performance Evaluation:**
- Calculate accuracy, precision, recall, F1-score, AUC
- Generate confusion matrices
- Per-class performance metrics

#### Results:
- Trained models with comprehensive performance metrics
- Model comparison analysis
- Best performing model identification

---

### Step 9: Cross-Validation and Overfitting Analysis

**Packages Used:** scikit-learn model_selection

#### Key Parameters:
- **Cross-validation:** StratifiedKFold, n_splits=5
- **Overfitting thresholds:** 0.05 (moderate), 0.10 (severe)
- **Learning curve training sizes:** 10 points from 10% to 100%

#### Processing Steps:

**1. Cross-Validation:**
- 5-fold stratified cross-validation
- Calculate training and validation scores

**2. Overfitting Assessment:**
- Analyze train-validation gaps
- Classify fitting quality (Good Fit, Moderate Overfitting, High Overfitting)

**3. Learning Curves:**
- Plot performance vs training set size
- Identify overfitting patterns

#### Results:
- Cross-validation performance metrics
- Overfitting analysis and recommendations
- Learning curve visualizations

---

## PHASE 3: Advanced Analysis and Visualization

### Step 10: Feature Importance Analysis

**Packages Used:** scikit-learn, pandas, numpy

#### Key Parameters:
- **Top features per model:** 20
- **Statistical significance threshold:** p < 0.05
- **Feature importance methods:** Gini importance, coefficient magnitudes, F-statistics

#### Processing Steps:

**1. Model-Specific Importance:**
- Extract feature weights from each trained model
- Tree-based: Gini importance
- Linear: Absolute coefficient values

**2. Cross-Model Analysis:**
- Identify consistently important features
- Rank features by discriminative power

#### Results:
- Comprehensive feature importance rankings
- Cross-model feature consistency analysis
- Top discriminative features identification

---

### Step 11: Cell Type-Specific Feature Analysis

**Packages Used:** scikit-learn feature_selection

#### Key Parameters:
- **Statistical test:** F-classification (ANOVA)
- **Significance threshold:** p < 0.05
- **Top features per cell type:** 15

#### Processing Steps:

**1. One-vs-All Analysis:**
- Binary classification for each cell type
- Statistical testing for discriminative features

**2. Cell Type Profiling:**
- Identify molecular signatures per cell type
- Rank features by statistical significance

#### Results:
- Cell type-specific molecular signatures
- Statistical significance scores for all features
- Discriminative feature profiles

---

### Step 12: Visualization Suite Creation

#### A. Feature Importance Heatmaps

**Packages Used:** matplotlib, seaborn, pandas

**Key Parameters:**
- Heatmap normalization: Column-wise (by model)
- Color scheme: viridis, plasma
- Feature display: Top 20 features

**Results:**
- Cross-model feature importance heatmap
- Cell type-specific feature heatmap
- RNA vs ATAC feature distribution analysis

---

#### B. Learning Curves and Performance Analysis

**Packages Used:** matplotlib, scikit-learn

**Key Parameters:**
- Training sizes: 8 points from 10% to 100%
- Cross-validation folds: 3
- Performance metric: Accuracy

**Results:**
- Learning curves for overfitting detection
- Training vs validation performance plots
- Model efficiency analysis

---

#### C. Network Analysis Visualization

**Packages Used:** networkx, matplotlib

**Key Parameters:**
- Feature correlation threshold: 0.5
- Network layout: Spring layout (k=2, iterations=50)
- Node sizes: Features=300, Cell types=500
- Top features included: 15

**Processing Steps:**

**1. Network Construction:**
- Features as nodes (potential regulatory elements)
- Cell types as nodes
- Correlation-based edges between features
- Association-based edges between features and cell types

**2. Visualization:**
- Color-coded nodes (features vs cell types)
- Edge styling (dashed for correlations, solid for associations)

**Results:**
- Feature-cell type association network
- Correlation-based feature relationships
- Network topology analysis

---

#### D. Performance Radar Charts

**Packages Used:** matplotlib (polar projection), numpy

**Parameters:**
- Metrics displayed: Accuracy, Precision, Recall, F1-Score, AUC
- Chart type: Polar coordinate system
- Color scheme: Distinct colors per model

**Results:**
- Multi-dimensional model comparison
- Performance trade-off visualization

---

#### E. Feature Distribution Analysis

**Packages Used:** seaborn, matplotlib, pandas

**Parameters:**
- Top features analyzed: 6 most important
- Plot type: Violin plots
- Grouping: By cell type

**Results:**
- Feature expression patterns across cell types
- Distribution shape analysis per cell type

---

### Step 13: Results Export and Documentation

**Packages Used:** json, pandas, openpyxl

#### Export Components:

**1. Model Artifacts:**
- Trained models (.pkl files)
- Label encoder and feature scaler
- Model parameters and configurations

**2. Analysis Results:**
- Performance metrics (CSV)
- Feature importance rankings (CSV)
- Cross-validation results (CSV)
- Cell type-specific analyses (CSV)

**3. Visualizations:**
- All plots in high-resolution PNG format
- Interactive visualizations (HTML)

**4. Comprehensive Report:**
- JSON format analysis summary
- Excel workbook with multiple sheets
- Markdown summary report

#### Results:
- Complete analysis documentation
- Reproducible result files
- Publication-ready visualizations

---

## Summary of Critical Parameters for Publication

### Data Processing Parameters
- **Quality Control:** min_genes=500, max_genes=7000, max_mito=20%
- **Normalization:** LogNormalize, scale_factor=10000
- **Feature Selection:** 2000 variable genes (RNA), 5000 peaks (ATAC)
- **Dimensionality:** PCA 1-15, UMAP default parameters

### Machine Learning Parameters
- **Train/Test Split:** 70/30, stratified by cell type
- **Class Balancing:** BorderlineSMOTE when max/min ratio > 2.0
- **Cross-Validation:** 5-fold StratifiedKFold
- **Overfitting Thresholds:** 0.05 (moderate), 0.10 (severe)

### Model-Specific Parameters
- **Random Forest:** 100 estimators, max_depth=10, balanced weights
- **XGBoost:** learning_rate=0.1, regularization (L1=0.1, L2=1.0)
- **SVM:** RBF kernel, C=1.0, balanced weights

### Statistical Analysis Parameters
- **Feature Selection:** F-classification, p<0.05
- **Network Analysis:** correlation threshold=0.5
- **Significance Testing:** ANOVA F-test for cell type discrimination

---

# Software Requirements 
##( you can check this READme.md file  --> for more details on how to setup your environment https://github.com/omicscodeathon/scatactf/blob/main/scripts/READme.md )

### R Environment (v4.5.1)
- Seurat v5.3.0
- Signac v1.15.0
- GenomicRanges
- EnsDb.Hsapiens.v75
- SingleR v2.0.0
- celldex
- hdf5r
- dplyr
- ggplot2
- biovizBase

### Python Environment (v3.10.17)
- pandas v2.2.3
- numpy v2.1.2
- scikit-learn v1.6.1
- xgboost v3.0.5
- imbalanced-learn v0.14.0
- matplotlib v3.10.1
- seaborn v0.13.2
- plotly v6.3.0
- networkx v3.3
- openpyxl v3.1.5
- joblib v1.4.2
- jupyter v1.1.1

---

## Contact

If you have any questions, you can contact team members .

---
