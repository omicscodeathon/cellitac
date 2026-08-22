# cellitac

**Cell type Identification using Transcription factor Analysis and Chromatin accessibility**

A pipeline for single-cell multiome (scATAC + scRNA) data that identifies cell types from **transcription-factor motif activity**. RNA is used only to derive labels; the classifier itself is trained purely on TF activity, so the model learns chromatin-level regulatory signal.

---

## Pipeline

| Stage | Steps | Tools |
|-------|-------|-------|
| **1. Preprocessing** (R) | multiome H5 → joint RNA+ATAC QC on shared barcodes | Seurat, Signac |
| | SingleR labels from the Monaco immune reference | SingleR, celldex |
| | JASPAR motif scan → chromVAR per-cell TF activity | JASPAR2020, motifmatchr, chromVAR |
| **2. Machine learning** (Python) | class composition, unsupervised feature cleaning | pandas, scikit-learn |
| | Logistic Regression, Random Forest, SVM, XGBoost | scikit-learn, xgboost, imbalanced-learn |
| | TF ↔ cell-type association (Mann-Whitney + BH-FDR + effect size) | scipy |
| | figures, tables and a JSON report | matplotlib, seaborn, networkx |

**Scope.** cellitac has been developed and tested on human **PBMC** multiome data, using the Monaco immune reference for SingleR labels. It works on any human tissue whose cell types are covered by that reference (blood, bone marrow, immune infiltrates).

**Genome builds:** hg38 (default) and hg19.

---

## Requirements

- Linux or macOS (Windows via WSL)
- Python 3.9 – 3.12
- Conda / Miniconda
- R ≥ 4.3 with the Bioconductor packages listed below (for the preprocessing stage). The ML stage runs on Python alone.

---

## Installation

Take the R packages from conda as pre-built binaries — do not let `BiocManager` compile them from source.

```bash
conda create -n cellitac -c conda-forge -c bioconda -y \
  python=3.11 rpy2 r-base=4.4 \
  r-seurat r-signac r-data.table \
  bioconductor-jaspar2020 bioconductor-tfbstools \
  bioconductor-motifmatchr bioconductor-chromvar \
  bioconductor-singler bioconductor-celldex \
  bioconductor-biovizbase bioconductor-rtracklayer \
  bioconductor-summarizedexperiment bioconductor-biocparallel \
  bioconductor-bsgenome.hsapiens.ucsc.hg38 \
  bioconductor-ensdb.hsapiens.v86

conda activate cellitac
pip install cellitac
