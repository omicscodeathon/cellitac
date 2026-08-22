# ==============================================================================
# cellitac — MERGED PREPROCESSING (single workflow, leakage-safe, TF-centric)
# ------------------------------------------------------------------------------
# Replaces the old R_Analysis_Part1.r (RNA) + R_Analysis_Part2.r (ATAC) split.
#
# WHAT THIS FIXES (directly addressing the ISMB-1088 reviews):
#   (A) LABEL LEAKAGE  -> RNA is used ONLY to derive cell-type labels (SingleR).
#                         RNA expression is NEVER written into the model feature
#                         matrix. The main model is trained on ATAC alone.
#   (B) NOT TF-CENTRIC -> Features are chromVAR per-cell TF *motif activity*
#                         (deviation z-scores) computed after JASPAR motif
#                         scanning of the peaks. No "top-5000 peaks" heuristic,
#                         no "assumed TF" features.
#   (C) GENOME BUILD   -> 10x PBMC 10k multiome is GRCh38/hg38. The old code
#                         used hg19/EnsDb.v75/blacklist_hg19, which would make
#                         motif matching wrong. This script uses hg38 throughout.
#                         >>> If YOUR fragments/peaks were aligned to a different
#                         build, change the GENOME BUILD block below to match. <<<
#   (D) BASELINES      -> RNA features are exported SEPARATELY (clearly labelled)
#                         so you can run the RNA-only and ATAC-only baselines the
#                         reviewers asked for. They are NEVER concatenated here.
#
# OUTPUT (all rows = the SAME matched cell barcodes, no cross-modality intersect):
#   cellitac_output/
#     |- cellitac_TF_activity.csv     <- X for the model: cells x TF motifs (chromVAR)
#     |- cell_labels.csv              <- y: cell_id, cell_type  (RNA-derived)
#     |- motif_to_TF_map.csv          <- JASPAR motif ID  ->  TF gene symbol
#     |- rna_baseline_features.csv    <- RNA expr (RNA-ONLY baseline; DO NOT merge)
#     |- atac_peak_baseline_lsi.csv   <- LSI embedding (ATAC peak baseline; optional)
#     |- multiome_processed.rds       <- full Seurat object (RNA + ATAC + chromvar)
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. LIBRARIES
# ------------------------------------------------------------------------------
# Install once (uncomment if needed):
   BiocManager::install(c("Signac","EnsDb.Hsapiens.v86","JASPAR2020","TFBSTools",
                          "BSgenome.Hsapiens.UCSC.hg38","motifmatchr","chromVAR",
                          "SingleR","celldex", "biovizBase"))
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)

  # ---- TF / motif stack (this is what makes the pipeline genuinely TF-centric)
  library(JASPAR2020)                    # TF motif database (PFMs)
  library(TFBSTools)                     # getMatrixSet()
  library(motifmatchr)                   # peak x motif matching (used by AddMotifs)
  library(chromVAR)                      # per-cell motif deviation scores
  library(BiocParallel)                  # backend control for chromVAR

  # ---- annotation (label side only)
  library(SingleR)
  library(celldex)
})

set.seed(42)

# ------------------------------------------------------------------------------
# PARALLEL BACKEND — FIX for the "sendMaster ... ignoring SIGPIPE" crash.
# RunChromVAR forks workers via BiocParallel; fork backends are unstable inside
# containers / Rscript / RStudio and the broken pipe kills the run. Force SERIAL
# execution = no forking, fully stable (a bit slower, fine for this cell count).
# ------------------------------------------------------------------------------
register(SerialParam())

# ------------------------------------------------------------------------------
# GENOME BUILD — set ONCE here and reuse. Must match how YOUR data was aligned.
# 10x "PBMC unsorted 10k" multiome demo data == GRCh38 (hg38).
# ------------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
GENOME_BSG  <- BSgenome.Hsapiens.UCSC.hg38
GENOME_NAME <- "hg38"
ENSDB       <- EnsDb.Hsapiens.v86
data("blacklist_hg38_unified", package = "Signac")   # hg38 ENCODE blacklist
BLACKLIST   <- blacklist_hg38_unified

getwd()
setwd("C:/Users/Niyonzima/Downloads/Cellitac/")
# ==============================================================================
# PATHS — everything (downloads AND results) lives under ONE folder, so nothing
# gets mixed up with the old 3k run.
# ==============================================================================
BASE_DIR <- "C:/Users/Niyonzima/Downloads/Cellitac/cellitac_fix"
DATA_DIR <- file.path(BASE_DIR, "data")      # raw 10x downloads
OUT_DIR  <- file.path(BASE_DIR, "results")   # pipeline outputs


# ==============================================================================
# STEP 0 — DOWNLOAD THE unsorted-10k DATASET (all files from ONE release)
# ------------------------------------------------------------------------------
# RELEASE: cell-arc / 2.0.0 / pbmc_unsorted_10k   (confirmed from the dataset page)
#
# The MERGED pipeline only reads 3 files, so we only download 3:
#   - filtered_feature_bc_matrix.h5  -> RNA + ATAC peak counts (peaks come from H5)
#   - atac_fragments.tsv.gz          -> QC (TSS/nucleosome) + chromVAR
#   - atac_fragments.tsv.gz.tbi      -> tabix index for the fragments (required)
# We deliberately DON'T grab atac_peaks.bed / atac_peak_annotation.tsv /
# per_barcode_metrics.csv: peaks are taken from the H5, annotations from EnsDb,
# and QC is computed from the data — so those files are never read here.
#
# All 3 come from one base URL => guaranteed same sample + same release.
# ==============================================================================
download_10x_unsorted_10k <- function(dest = DATA_DIR, overwrite = FALSE) {
  options(timeout = 7200)   # fragments file is large (GBs) -> long timeout
  if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)

  base   <- "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k"
  sample <- "pbmc_unsorted_10k"

  # ONLY the files the pipeline actually consumes, ALL from release 2.0.0:
  files <- c(
    "filtered_feature_bc_matrix.h5",
    "atac_fragments.tsv.gz",
    "atac_fragments.tsv.gz.tbi"
  )

  cat("[0] Downloading pbmc_unsorted_10k (cell-arc/2.0.0) -> ", dest, "\n")
  cat("    source release:", base, "\n")

  for (f in files) {
    url  <- file.path(base, paste0(sample, "_", f))
    out  <- file.path(dest, paste0(sample, "_", f))
    if (file.exists(out) && !overwrite && file.size(out) > 0) {
      cat("    [skip] already present:", basename(out),
          sprintf("(%.1f MB)\n", file.size(out) / 1e6))
      next
    }
    cat("    [get ] ", basename(out), " ...\n")
    # mode='wb' is mandatory for binary files (.h5/.gz) on all platforms
    ok <- tryCatch(
      { download.file(url, destfile = out, mode = "wb", quiet = TRUE); TRUE },
      error = function(e) { cat("       FAILED:", conditionMessage(e), "\n"); FALSE }
    )
    if (!ok || !file.exists(out) || file.size(out) == 0)
      stop("Download failed or empty file: ", url,
           "\n  -> check the URL/version on the dataset page if 10x moved it.")
  }

  # --- VERIFY: every expected file exists, is non-empty, same prefix ----------
  cat("[0] Verifying all files are from the SAME sample/release...\n")
  report <- data.frame(
    file = paste0(sample, "_", files),
    size_MB = round(file.size(file.path(dest, paste0(sample, "_", files))) / 1e6, 2),
    exists  = file.exists(file.path(dest, paste0(sample, "_", files)))
  )
  print(report)
  if (!all(report$exists) || any(report$size_MB == 0))
    stop("Some files are missing or empty — do NOT proceed until fixed.")
  # all share the same sample prefix => guaranteed same release/batch
  stopifnot(all(grepl(paste0("^", sample, "_"), report$file)))
  cat("    OK:", nrow(report), "/", nrow(report),
      "files present, same release (", sample, ", cell-arc/2.0.0 ).\n")

  invisible(dest)
}


# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================
run_cellitac_preprocessing <- function(
    input_dir  = DATA_DIR,
    output_dir = OUT_DIR,
    h5_name        = "pbmc_unsorted_10k_filtered_feature_bc_matrix.h5",
    fragments_name = "pbmc_unsorted_10k_atac_fragments.tsv.gz"
) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  h5_file        <- file.path(input_dir, h5_name)
  fragments_file <- file.path(input_dir, fragments_name)
  stopifnot(file.exists(h5_file), file.exists(fragments_file))

  # ----------------------------------------------------------------------------
  # 1. LOAD MULTIOME H5 ONCE  (gives BOTH modalities on the SAME barcodes)
  #    -> this is why we merge the two old scripts: no more intersect() that
  #       silently drops cells between RNA and ATAC objects.
  # ----------------------------------------------------------------------------
  cat("[1] Loading multiome H5 (RNA + ATAC together)...\n")
  data_10x   <- Read10X_h5(h5_file)
  rna_counts <- data_10x[["Gene Expression"]]
  atac_counts<- data_10x[["Peaks"]]

  # Build RNA assay first
  obj <- CreateSeuratObject(counts = rna_counts, assay = "RNA", project = "cellitac")

  # ----------------------------------------------------------------------------
  # 2. ATTACH ATAC AS A ChromatinAssay  (peaks restricted to standard chromosomes)
  # ----------------------------------------------------------------------------
  cat("[2] Building ChromatinAssay on", GENOME_NAME, "...\n")
  grange_counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  std_idx       <- as.vector(seqnames(grange_counts) %in% standardChromosomes(grange_counts))
  atac_counts   <- atac_counts[std_idx, ]

  annotations <- GetGRangesFromEnsDb(ensdb = ENSDB)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- GENOME_NAME

  chrom_assay <- CreateChromatinAssay(
    counts    = atac_counts,
    sep       = c(":", "-"),
    fragments = fragments_file,
    genome    = GENOME_NAME,
    annotation= annotations,
    min.cells = 10
  )
  # keep only cells present in BOTH modalities
  common <- intersect(colnames(obj), colnames(chrom_assay))
  obj <- subset(obj, cells = common)
  obj[["ATAC"]] <- subset(chrom_assay, cells = common)

  cat("    cells:", ncol(obj),
      "| RNA genes:", nrow(obj[["RNA"]]),
      "| ATAC peaks:", nrow(obj[["ATAC"]]), "\n")

  # ----------------------------------------------------------------------------
  # 3. JOINT QC  (a cell must pass BOTH RNA and ATAC thresholds)
  # ----------------------------------------------------------------------------
  cat("[3] Joint QC (RNA + ATAC)...\n")
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  DefaultAssay(obj) <- "ATAC"
  obj <- NucleosomeSignal(obj)
  # fast = TRUE: computes the TSS score WITHOUT storing the full per-position
  # matrix. Much lower memory + far faster on 12k cells / 2.9GB fragments, and
  # the score is all we need for QC filtering (fast=FALSE is only needed to draw
  # the TSS enrichment plot, which we don't use here).
  obj <- TSSEnrichment(obj, fast = TRUE)
  obj$blacklist_ratio <- FractionCountsInRegion(obj, assay = "ATAC", regions = BLACKLIST)
  # NOTE: pct_reads_in_peaks was removed — it re-read the 2.9GB fragments file
  # TWICE and was never used in the filter, so it only wasted time/memory.

  n0 <- ncol(obj)

  # --- DIAGNOSTIC: how many cells FAIL each threshold individually -------------
  # (tells you which single filter is responsible if the drop is too large)
  md <- obj@meta.data
  fail <- c(
    nFeature_RNA_lo   = sum(md$nFeature_RNA  <= 500,   na.rm = TRUE),
    nFeature_RNA_hi   = sum(md$nFeature_RNA  >= 7000,  na.rm = TRUE),
    nCount_RNA_lo     = sum(md$nCount_RNA    <= 1000,  na.rm = TRUE),
    percent_mt_hi     = sum(md$percent.mt    >= 20,    na.rm = TRUE),
    nCount_ATAC_lo    = sum(md$nCount_ATAC   <= 1000,  na.rm = TRUE),
    nCount_ATAC_hi    = sum(md$nCount_ATAC   >= 100000,na.rm = TRUE),
    TSS_enrich_lo     = sum(md$TSS.enrichment<= 2,     na.rm = TRUE),
    nucleosome_hi     = sum(md$nucleosome_signal >= 4, na.rm = TRUE),
    blacklist_hi      = sum(md$blacklist_ratio   >= 0.05, na.rm = TRUE)
  )
  cat("    cells failing EACH threshold (out of", n0, "):\n")
  print(fail)

  obj <- subset(
    obj,
    subset = nFeature_RNA  > 500   & nFeature_RNA  < 7000 &
             nCount_RNA     > 1000  & percent.mt    < 20    &
             nCount_ATAC    > 1000  & nCount_ATAC   < 100000 &
             TSS.enrichment > 2     & nucleosome_signal < 4  &
             blacklist_ratio < 0.05
  )
  cat("    cells:", n0, "->", ncol(obj),
      "(", round(ncol(obj)/n0*100, 1), "% kept )\n")

  # ----------------------------------------------------------------------------
  # 4. RNA SIDE  ==>  LABELS ONLY  (this side never becomes a model feature)
  # ----------------------------------------------------------------------------
  cat("[4] RNA processing -> SingleR labels (LABELS ONLY, not features)...\n")
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj) |>
         FindVariableFeatures(nfeatures = 2000) |>
         ScaleData() |>
         RunPCA(verbose = FALSE) |>
         FindNeighbors(dims = 1:15) |>
         FindClusters(resolution = 0.5, verbose = FALSE) |>
         RunUMAP(dims = 1:15, verbose = FALSE)

  # Monaco immune reference gives clean PBMC immune labels. HPCA is a general
  # atlas that injects bone-marrow progenitor labels (CMP/HSC/Pre-B_CD34-) that
  # don't belong in blood -> those spurious 2-cell classes disappear with Monaco.
  ref   <- celldex::MonacoImmuneData()
  sce   <- as.SingleCellExperiment(obj)
  preds <- SingleR(sce, ref = ref, labels = ref$label.main)
  obj$cell_type <- preds$labels[match(colnames(obj), rownames(preds))]
  cat("    label distribution:\n"); print(table(obj$cell_type))

  # ----------------------------------------------------------------------------
  # 5. ATAC SIDE  ==>  TF MOTIF ACTIVITY (the real TF-centric features)
  #    JASPAR motif scan -> chromVAR per-cell deviation z-scores.
  # ----------------------------------------------------------------------------
  cat("[5] ATAC processing -> JASPAR scan + chromVAR TF activity...\n")
  DefaultAssay(obj) <- "ATAC"
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = "q0")
  obj <- RunSVD(obj)                                  # LSI (used for ATAC baseline)
  obj <- RunUMAP(obj, reduction = "lsi", dims = 2:30, verbose = FALSE)

  # 5a. Pull JASPAR2020 CORE vertebrate motifs for Homo sapiens (taxid 9606)
  pfm <- getMatrixSet(JASPAR2020, opts = list(
    collection = "CORE", species = 9606, all_versions = FALSE))

  # 5b. Scan every peak for every motif (peak x motif binary match matrix)
  obj <- AddMotifs(obj, genome = GENOME_BSG, pfm = pfm)

  # 5c. chromVAR: per-cell, per-motif deviation z-score = quantitative TF activity
  #obj <- RunChromVAR(obj, genome = GENOME_BSG)
  # --- REPLACE THE LINE BELOW WITH THIS BLOCK ---
  library(chromVAR)
  library(SummarizedExperiment)
  
  # Prepare the SummarizedExperiment
  se <- SummarizedExperiment(
    assays = list(counts = GetAssayData(obj, assay = "ATAC", layer = "counts")),
    rowRanges = granges(obj[["ATAC"]]),
    colData = obj@meta.data
  )
  
  # Calculate GC bias and then the deviations
  se <- addGCBias(se, genome = GENOME_BSG)
  motif_ix <- obj[["ATAC"]]@motifs@data
  deviations_obj <- computeDeviations(object = se, annotations = motif_ix)
  
  # Add to object as "chromvar" assay
  obj[["chromvar"]] <- CreateAssayObject(counts = assay(deviations_obj, "z"))
  DefaultAssay(obj) <- "chromvar"
  
  # ----------------------------------------------
  
  # ----------------------------------------------------------------------------
  # 6. EXPORT  (everything aligned to the SAME final barcodes)
  # ----------------------------------------------------------------------------
  cat("[6] Exporting matched matrices...\n")
  cells <- colnames(obj)

  ## 6a. MODEL FEATURES (X): chromVAR TF activity, cells x motifs
  DefaultAssay(obj) <- "chromvar"
  tf_activity <- t(as.matrix(GetAssayData(obj, layer = "data")))   # cells x motifs
  tf_activity <- tf_activity[cells, ]

  ## 6b. Map JASPAR motif IDs -> readable TF gene symbols (interpretability)
  motif_obj <- Motifs(obj[["ATAC"]])
  motif_map <- data.frame(
    motif_id = colnames(tf_activity),
    TF       = unlist(motif_obj@motif.names[colnames(tf_activity)]),
    row.names = NULL
  )
  # Use TF symbols as the exported column names so feature importances are readable
  colnames(tf_activity) <- make.unique(motif_map$TF)

  write.csv(tf_activity, file.path(output_dir, "cellitac_TF_activity.csv"))
  write.csv(motif_map,   file.path(output_dir, "motif_to_TF_map.csv"), row.names = FALSE)

  ## 6c. LABELS (y) — RNA-derived, matched to the SAME cells
  labels <- data.frame(cell_id = cells, cell_type = obj$cell_type[cells])
  labels <- labels[!is.na(labels$cell_type), ]
  write.csv(labels, file.path(output_dir, "cell_labels.csv"), row.names = FALSE)

  ## 6d. RNA-ONLY BASELINE features — exported SEPARATELY, never concatenated.
  ##     (Lets you run the RNA-only baseline the reviewers requested.)
  DefaultAssay(obj) <- "RNA"
  vg  <- head(VariableFeatures(obj), 2000)
  rna <- t(as.matrix(GetAssayData(obj, layer = "data")[vg, ]))[cells, ]
  write.csv(rna, file.path(output_dir, "rna_baseline_features.csv"))

  ## 6e. ATAC-PEAK BASELINE — LSI embedding (label-free peak baseline, optional)
  lsi <- Embeddings(obj, "lsi")[cells, 2:30]
  write.csv(lsi, file.path(output_dir, "atac_peak_baseline_lsi.csv"))

  saveRDS(obj, file.path(output_dir, "multiome_processed.rds"))

  cat("\n=== DONE ===\n")
  cat("Model features (ATAC TF activity):", nrow(tf_activity), "cells x",
      ncol(tf_activity), "TF motifs\n")
  cat("Files written to:", normalizePath(output_dir), "\n")
  cat("REMINDER: train the main model on cellitac_TF_activity.csv ONLY.\n")
  cat("          rna_baseline_features.csv is for the RNA-only baseline — do NOT merge.\n")

  invisible(obj)
}

# ==============================================================================
# RUN  (download first, then process — all under /home/rhamed/cellitac_fix)
# ==============================================================================
download_10x_unsorted_10k()                 # Step 0: fetch + verify same release
obj <- run_cellitac_preprocessing()         # Steps 1-6: process the unsorted 10k

