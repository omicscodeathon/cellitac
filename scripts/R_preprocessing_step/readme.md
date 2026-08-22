Run this before the Python ML pipeline.

Downloads the 10x PBMC unsorted-10k multiome set, does joint RNA+ATAC QC,
pulls cell-type labels from SingleR (Monaco reference), then computes
chromVAR TF motif activity per cell (JASPAR, hg38 by default).

Writes to `cellitac_output/`:
- `cellitac_TF_activity.csv` — TF activity per cell, this is the model input
- `cell_labels.csv` — cell types, this is the target
- `motif_to_TF_map.csv` — JASPAR motif ID to TF gene symbol

`cellitac_ML_pipeline_v3.py` won't run without these three files.
