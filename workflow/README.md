# cellitac — Analysis Workflow

Cell type identification from single-cell multiome (scATAC + scRNA) data using
**transcription-factor motif activity** as the sole feature space. RNA is used
only to derive cell-type labels; it never enters the classifier.

For installation and CLI usage see the [README](./README.md).

---

## Design principles

1. **No label leakage.** RNA gene expression is used to compute SingleR labels
   and is then dropped. The classifier sees only ATAC-derived TF motif activity,
   so labels and features come from independent modalities.
2. **TF-centric features.** Peaks are scanned against JASPAR2020 CORE
   vertebrate motifs and summarised per cell with chromVAR deviation z-scores,
   so features are directly interpretable as putative TF activity.
3. **All supervised steps are fitted in-fold.** Scaling, feature selection and
   (optionally) resampling live inside a scikit-learn `Pipeline`, so nothing
   fitted on training folds ever sees the test set.
4. **Reproducibility.** A single random seed (42) is threaded through every
   stochastic step.

---

## Stage 1 — Preprocessing (R via rpy2)

Runs from a single `.h5` + fragments file. Genome build is configurable
(`--genome hg38` default, or `hg19`).

**Packages:** Seurat, Signac, celldex, SingleR, JASPAR2020, TFBSTools,
motifmatchr, chromVAR, BiocParallel.

### 1a. Load and build a joint multiome object

- The 10x `.h5` is loaded once, giving RNA counts and Peak counts on the
  **same barcodes** — no cross-modality intersection is needed.
- Signac `ChromatinAssay` is created with `min.cells = 10`, fragments file
  attached for downstream QC.

### 1b. Joint QC

A cell is kept only if it passes both RNA and ATAC thresholds simultaneously:

| Metric | Threshold |
|---|---|
| `nFeature_RNA` | > 500 and < 7000 |
| `nCount_RNA` | > 1000 |
| `percent.mt` | < 20% |
| `nCount_ATAC` | > 1000 and < 100 000 |
| `TSS.enrichment` | > 2 |
| `nucleosome_signal` | < 4 |

### 1c. RNA processing → labels only

Standard Seurat: `NormalizeData` → `FindVariableFeatures` (2000 genes) →
`ScaleData` → `RunPCA`. Labels are assigned with **SingleR** against the
**Monaco immune reference** (`celldex::MonacoImmuneData()`, `label.main`).

RNA outputs are used only to produce `cell_labels.csv`. They are not exported
as features.

### 1d. ATAC processing → TF motif activity

- `RunTFIDF` → `FindTopFeatures(min.cutoff = "q0")` → `RunSVD` → LSI dims 2–30
  for UMAP.
- JASPAR2020 CORE vertebrate motifs for *Homo sapiens*
  (`getMatrixSet(JASPAR2020, opts = list(collection = "CORE", species = 9606))`).
- `matchMotifs` (motifmatchr) against the peak set using the chosen
  `BSgenome.Hsapiens.UCSC.*`.
- `chromVAR::computeDeviations` on a `SerialParam()` backend gives the
  **per-cell per-motif deviation z-score** — the model's feature space.

### Stage 1 outputs

```
preprocessing/
├── cellitac_TF_activity.csv    # X — cells × TF motifs (chromVAR z-scores)
├── cell_labels.csv             # y — cell_id, cell_type
├── motif_to_TF_map.csv         # JASPAR motif ID → TF gene symbol
└── multiome_processed.rds      # full Seurat object (checkpoint)
```

---

## Stage 2 — Machine learning (Python)

Runs on the CSVs alone; needs no R. Implemented as `CellitacPipeline` in
`cellitac.mainModel`.

**Packages:** scikit-learn, xgboost, imbalanced-learn, scipy, matplotlib,
seaborn, networkx.

### 2a. Load and guard against leakage

- `cellitac_TF_activity.csv` and `cell_labels.csv` are loaded and intersected
  on barcodes.
- **Leakage tripwire:** if any column name matches a known RNA marker
  (`CD3E`, `MS4A1`, `LYZ`, `NKG7`, `CD14`, `CD8A`, `MS4A1`, `PPBP`, …) the
  run aborts with a clear error.
- Classes with fewer than `MIN_CELLS = 50` cells are dropped so cross-
  validation folds always contain the class.

### 2b. Unsupervised feature cleaning (no labels used)

- Zero-variance motifs removed.
- Motifs with absolute Spearman correlation > `CORR_CUTOFF = 0.95` collapsed
  to one representative.

### 2c. Stratified train/test split — **before any supervised step**

`train_test_split(test_size = 0.2, stratify = y, random_state = 42)`.

### 2d. Per-model pipeline

Every model is wrapped in an `imblearn.Pipeline`:

```
StandardScaler → SelectKBest(f_classif, k = min(MAX_FEATURES, n_motifs))
              → [optional SMOTE — off by default]
              → classifier
```

`MAX_FEATURES = 150`. Because selection lives inside the pipeline, the
K best motifs are recomputed on the training fold of every CV split — no
leakage from selection.

### 2e. Models

| Model | Key hyperparameters |
|---|---|
| Logistic Regression | L2, C = 0.1, `class_weight = "balanced"`, `max_iter = 3000` |
| Random Forest | 300 trees, `max_depth = 6`, `min_samples_leaf = 15`, `max_features = "sqrt"`, `class_weight = "balanced"` |
| XGBoost | 300 trees, `max_depth = 3`, `lr = 0.05`, `subsample = 0.8`, `colsample_bytree = 0.8`, `reg_alpha = 1.0`, `reg_lambda = 5.0` |
| SVM | RBF, C = 0.5, `gamma = "scale"`, `class_weight = "balanced"`, `probability = True` |

Class imbalance is handled with `class_weight = "balanced"` by default; SMOTE
is available (`USE_SMOTE = True`) but not on by default because synthetic
cells in a leakage-safe setup add risk with little gain.

### 2f. Evaluation

- **5-fold stratified CV** (`StratifiedKFold(n_splits = 5, shuffle = True,
  random_state = 42)`) on the training set → train-vs-validation gap → an
  overfitting verdict per model.
- **Held-out test set metrics:** balanced accuracy, macro F1, weighted AUC
  (one-vs-rest for multiclass).
- **Learning curves** on 6 training sizes with 3-fold CV.

### 2g. TF ↔ cell-type association

For every TF motif and every class:

- One-vs-rest **Mann-Whitney U test** on chromVAR z-scores (no normality
  assumed).
- **Benjamini-Hochberg FDR** across all TF × class tests (`FDR_ALPHA = 0.05`).
- **Effect size:** Cohen's *d* (`≥ 0.50`) and rank-biserial AUC (`≥ 0.60`).
- A TF is called significant for a class only if all three pass.

### 2h. Interpretation outputs

- Per-model top-20 TFs (native importance for tree/linear models, permutation
  importance for kernel SVM).
- **Cell ↔ TF network** (`fig07`, `table08`): edges only where the TF is
  significant for the class after BH-FDR + effect-size filters. TF-TF edges
  disabled by default.

### Stage 2 outputs

```
ml_results/
├── cellitac_ml_report.json
├── fig01_umap_TFactivity_before_training.png
├── fig02_class_composition_pies.png
├── fig03_class_imbalance_handling.png
├── fig04_model_comparison.png
├── fig05_confusion_matrices.png
├── fig06_importance_<model>.png
├── fig07_TF_network_<model>.png
├── fig08_learning_curves.png
├── fig09_umap_tsne_after_training.png
├── table01_class_composition.csv
├── table02_model_performance.csv
├── table03_accuracy_and_overfitting.csv
├── table04_per_class_metrics.csv
├── table05_per_class_recall.csv
├── table06_top20_TF_all_models.csv
├── table07_TF_celltype_associations_full.csv
├── table08_network_edges_<model>.csv
└── table09_learning_curve.csv
```

---

## Scope and limitations

- Tested on human **PBMC** multiome data. Any human tissue whose cell types
  are covered by the Monaco immune reference should work (blood, bone marrow,
  immune infiltrates). Other tissues require swapping the SingleR reference
  in the R stage.
- Human only: JASPAR is queried with `species = 9606` and genome resources
  are `BSgenome.Hsapiens` and `EnsDb.Hsapiens`.
- Genome build must match how the fragments file was aligned; hg38 default,
  hg19 supported.

---

## Contact

Open an issue on the [GitHub repository](https://github.com/omicscodeathon/cellitac).
