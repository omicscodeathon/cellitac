# ================================================================================
# Part two : ATAC DATA PROCESSING & FEATURE ENGINEERING  
# ================================================================================

# Installing packages
#BiocManager::install("Signac")
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("biovizBase")
library(Signac)
library(GenomicRanges)
library(data.table)
library(EnsDb.Hsapiens.v75) # Use appropriate annotation package for your species
library(biovizBase)

atac_processing_team2 <- function() {
  
  cat("Part two 2: ATAC Data Processing Started ===\n")
  
  # Setup the directories
  input_dir1 <- "~/singlecell/ATAC/team1_rna_output/"
  input_dir2 <- "~/singlecell/ATAC/"
  output_dir_team2 <- "team2_atac_output"
  if (!dir.exists(output_dir_team2)) dir.create(output_dir_team2, recursive = TRUE)
  
  # ----------------------------------------
  # LOAD ATAC DATA FILES...from same URL
  # ----------------------------------------
  
  cat("Loading ATAC-specific files...\n")
  
  # Main files for ATAC processing
  rds_file <- file.path(input_dir1, "pbmc_rna_processed.rds")
  fragments_file <- file.path(input_dir2, "pbmc_unsorted_10k_atac_fragments.tsv.gz")  
  peaks_file <- file.path(input_dir2, "pbmc_unsorted_10k_atac_peaks.bed")            
  
  # Load rds data and extract ATAC portion
  data_10x <- readRDS(rds_file)
  
  # Loading the peaks file (BED format)
  peaks <- read.table(
    file = peaks_file,
    col.names = c("chr", "start", "end")
  )
  
  # Convert to GRanges object
  granges <- makeGRangesFromDataFrame(peaks)
  
  # Count fragments per cell per peak
  atac_counts <- FeatureMatrix(
    fragments = fragments,
    features = granges, 
    cells = NULL # Will determine cells automatically
  )
  
  #cat("ATAC data created with dimensions:", dim(atac_counts), "(peaks x cells)\n")
  cat("ATAC data dimensions:", nrow(atac_counts), "peaks x", ncol(atac_counts), "cells\n")

  # Get annotations for the appropriate genome
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
  seqlevelsStyle(annotations) <- 'UCSC'
  
  
  #if (is.list(data_10x)) {
  #  atac_counts <- data_10x$Peaks
  #  cat("Extracted ATAC portion from multiome data\n")
  #} else {
  #  stop("Expected multiome data with separate RNA and ATAC components")
  #}
  
  
  # ----------------------------------------
  # CREATE CHROMATIN ASSAY...
  # ----------------------------------------
  
  cat("Creating ChromatinAssay object...\n")
  
  # Create ChromatinAssay
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    annotation = annotations,
    fragments = fragments_file,
    genome = "hg19",
    min.cells = 10,
    min.features = 200
  )
  
  # Create Seurat object with ChromatinAssay... 
  pbmc_atac <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    project = "PBMC_ATAC"
  )
  
  cat("ChromatinAssay created with", ncol(pbmc_atac), "cells and", nrow(pbmc_atac), "peaks\n")
  
  # ----------------------------------------
  # ATAC QC METRICS CALCULATION
  # ----------------------------------------
  
  cat("Calculating ATAC-specific QC metrics...\n")
  
  # Add blacklist information
  blacklist <- get(data(list = "blacklist_hg19"))
  pbmc_atac$blacklist_fraction <- FractionCountsInRegion(
    object = pbmc_atac, 
    regions = blacklist
  )
  
  
  
  # Create the fragment object
  fragments <- CreateFragmentObject(
    path = fragments_file
  )
  
  # 1. Calculate total fragments per cell (passed_filters)
  cat("Calculating total fragments per cell...\n")
  total_fragments <- CountFragments(
    fragments = fragments_file,  # Use the file path directly
    cells = colnames(pbmc_atac)
  )
  
  # Add to metadata
  pbmc_atac$passed_filters <- total_fragments$frequency_count[match(colnames(pbmc_atac), 
                                                                    total_fragments$CB)]
  
  # 2. Calculate fragments in peaks
  cat("Calculating fragments in peaks...\n")
  # Get the genomic ranges of peaks from your assay
  peak_ranges <- granges(pbmc_atac)
  
  
  # Count fragments in peaks
  fragments_in_peaks <- CountFragments(
    fragments = fragments_file,
    cells = colnames(pbmc_atac)
  )
  
  pbmc_atac$peak_region_fragments <- fragments_in_peaks$frequency_count[match(colnames(pbmc_atac), fragments_in_peaks$CB)]
  
  # 
  cat("Calculating derived QC metrics...\n")
  pbmc_atac$pct_reads_in_peaks <- pbmc_atac$peak_region_fragments / pbmc_atac$passed_filters * 100
  pbmc_atac$blacklist_ratio <- pbmc_atac$blacklist_fraction / pbmc_atac$peak_region_fragments
  
  # ----------------------------------------
  # VERIFY THE METRICS
  # ----------------------------------------
  cat("QC metrics added to object:\n")
  print(head(pbmc_atac@meta.data))
  
  # Check which columns are now available
  cat("Available metadata columns:\n")
  print(colnames(pbmc_atac@meta.data))
  
  # Summary statistics
  summary(pbmc_atac@meta.data[, c("peak_region_fragments", "pct_reads_in_peaks", 
                                  "blacklist_ratio", "nucleosome_signal", "TSS.enrichment")])
  
  # See how many cells would pass typical filters
  cat("Cells passing typical thresholds:\n")
  cat("peak_region_fragments > 1000:", sum(pbmc_atac$peak_region_fragments > 1000), "\n")
  cat("pct_reads_in_peaks > 15:", sum(pbmc_atac$pct_reads_in_peaks > 15), "\n")
  cat("blacklist_ratio < 0.05:", sum(pbmc_atac$blacklist_ratio < 0.05), "\n")
  cat("nucleosome_signal < 4:", sum(pbmc_atac$nucleosome_signal < 4), "\n")
  cat("TSS.enrichment > 2:", sum(pbmc_atac$TSS.enrichment > 2), "\n")
  
  
  # Calculate QC metrics
  pbmc_atac <- NucleosomeSignal(object = pbmc_atac)
  pbmc_atac <- TSSEnrichment(object = pbmc_atac, fast = FALSE)
  #colnames(pbmc_atac@meta.data)
 
  # ----------------------------------------
  # ATAC QUALITY FILTERING.....
  # ----------------------------------------
  
  cat("Applying ATAC quality filters...\n")
  
  n_cells_before <- ncol(pbmc_atac)
  n_peaks_before <- nrow(pbmc_atac)
  
  # Define ATAC-specific thresholds
  min_peak_count <- 1000     # Min peak counts per cell
  max_peak_count <- 100000   # Max peak counts (remove potential doublets)
  min_tss_enrichment <- 2    # Min TSS enrichment
  max_nucleosome_signal <- 4 # Max nucleosome signal
  min_pct_in_peaks <- 15     # Min percentage of reads in peaks
  max_blacklist_ratio <- 0.05
  
  # Apply filters
  pbmc_atac <- subset(
    x = pbmc_atac,
    subset = peak_region_fragments > min_peak_count &
             peak_region_fragments < max_peak_count &
             pct_reads_in_peaks > min_pct_in_peaks &
             blacklist_ratio < max_blacklist_ratio &
             nucleosome_signal < max_nucleosome_signal &
             TSS.enrichment > min_tss_enrichment
  )
  
  cat("After ATAC QC filtering:\n")
  cat("  - Cells:", n_cells_before, "->", ncol(pbmc_atac), 
      "(", round(ncol(pbmc_atac)/n_cells_before*100, 1), "%)\n")
  cat("  - Peaks:", n_peaks_before, "->", nrow(pbmc_atac), 
      "(", round(nrow(pbmc_atac)/n_peaks_before*100, 1), "%)\n")
  
  # ----------------------------------------
  # ATAC DATA NORMALIZATION.....
  # ----------------------------------------
  
  cat("Normalizing ATAC data...\n")
  
  # TF-IDF normalization
  pbmc_atac <- RunTFIDF(pbmc_atac)
  
  # Feature selection - find top variable peaks
  pbmc_atac <- FindTopFeatures(pbmc_atac, min.cutoff = 'q0')
  
  # ----------------------------------------
  # DIMENSIONALITY REDUCTION FOR ATAC....UMAP LSI
  # ----------------------------------------
  
  cat("Performing ATAC dimensionality reduction...\n")
  
  # LSI (Latent Semantic Indexing) for ATAC
  pbmc_atac <- RunSVD(pbmc_atac)
  
  # Non-linear dimensionality reduction
  pbmc_atac <- RunUMAP(object = pbmc_atac, reduction = 'lsi', dims = 2:30)
  
  # Find neighbors and clusters
  pbmc_atac <- FindNeighbors(object = pbmc_atac, reduction = 'lsi', dims = 2:30)
  pbmc_atac <- FindClusters(object = pbmc_atac, verbose = FALSE, algorithm = 3)
  
  # ----------------------------------------
  # ATAC FEATURE PREPARATION....
  # ----------------------------------------
  
  cat("Preparing ATAC features for ML...\n")
  
  # Get top 5000 most accessible peaks
  peak_accessibility <- Matrix::rowSums(pbmc_atac@assays$ATAC@counts > 0)
  top_peaks <- names(sort(peak_accessibility, decreasing = TRUE))[1:5000]
  
  # Get normalized ATAC data for these peaks
  atac_data_normalized <- GetAssayData(pbmc_atac, assay = "ATAC", slot = "data")
  atac_features <- as.data.frame(as.matrix(t(atac_data_normalized[top_peaks, ])))
  
  # ----------------------------------------
  # EXPORT STEP 2 RESULTS...
  # ----------------------------------------
  
  cat("Exporting STEP 2 TEAM 2 results...\n")
  
  # Export ATAC feature matrix
  write.csv(atac_features, file.path(output_dir_team2, "atac_features_5000.csv"))
  
  # Export ATAC metadata
  write.csv(pbmc_atac@meta.data, file.path(output_dir_team2, "atac_metadata.csv"))
  
  # Export peak informationnn..
  peak_info <- data.frame(
    peak = top_peaks,
    accessibility_score = peak_accessibility[top_peaks]
  )
  write.csv(peak_info, file.path(output_dir_team2, "top_peaks_info.csv"), row.names = FALSE)
  
  # Create ATAC QC plots
  p1 <- VlnPlot(
    object = pbmc_atac,
    features = c('pct_reads_in_peaks', 'peak_region_fragments', 
                 'TSS.enrichment', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 4
  )
  
  p2 <- DimPlot(pbmc_atac, reduction = "umap", label = TRUE)
  
  ggsave(file.path(output_dir_team2, "atac_qc_plots.png"), p1, width = 16, height = 6)
  ggsave(file.path(output_dir_team2, "atac_clusters_umap.png"), p2, width = 10, height = 8)
  
  # Save ATAC Seurat object
  saveRDS(pbmc_atac, file.path(output_dir_team2, "pbmc_atac_processed.rds"))
  
  cat("step team 2 ATAC PROCESSING COMPLETED ===\n")
  cat("Output files saved in:", output_dir_team2, "\n")
  
  return(pbmc_atac)
}

# ================================================================================
# FINAL INTEGRATION: COMBINE TEAM 1 & 2 OUTPUTS FOR PYTHON ...for ML model building
# ================================================================================

combine_for_python <- function() {
  
  cat(" COMBINING TEAM OUTPUTS FOR PYTHON ML ..\n")
  
  # Setup directories
  output_dir_python <- "python_ready_data"
  if (!dir.exists(output_dir_python)) dir.create(output_dir_python, recursive = TRUE)
  
  # ----------------------------------------
  # LOAD TEAMS OUTPUTS ...
  # ----------------------------------------
  
  # Load RNA features from Team 1
  rna_features <- read.csv("team1_rna_output/rna_features_2000.csv", row.names = 1)
  rna_metadata <- read.csv("team1_rna_output/rna_metadata.csv", row.names = 1)
  cell_labels <- read.csv("team1_rna_output/cell_type_labels.csv")
  
  # Load ATAC features from Team 2  
  atac_features <- read.csv("team2_atac_output/atac_features_5000.csv", row.names = 1)
  atac_metadata <- read.csv("team2_atac_output/atac_metadata.csv", row.names = 1)
  
  # -------------------------------------------
  # MATCH CELLS BETWEEN MODALITIES...
  # ----------------------------------------
  
  # Find common cells between RNA and ATAC
  common_cells <- intersect(rownames(rna_features), rownames(atac_features))
  
  cat("Common cells between RNA and ATAC:", length(common_cells), "\n")
  
  # Subset to common cells
  rna_features_matched <- rna_features[common_cells, ]
  atac_features_matched <- atac_features[common_cells, ]
  
  # ----------------------------------------
  # CREATE COMBINED FEATURE MATRIX
  # ----------------------------------------
  
  # Combine RNA and ATAC features
  combined_features <- cbind(rna_features_matched, atac_features_matched)
  
  cat("Combined feature matrix dimensions:", nrow(combined_features), "cells x", 
      ncol(combined_features), "features\n")
  
  # ----------------------------------------
  # PREPARE LABELS AND METADATA
  # ----------------------------------------
  
  # match cell labels to common cells
  cell_labels_matched <- cell_labels[cell_labels$cell_id %in% common_cells, ]
  
  # Create final labels dataframe
  final_labels <- data.frame(
    cell_id = common_cells,
    cell_type = cell_labels_matched$cell_type[match(common_cells, cell_labels_matched$cell_id)]
  )
  
  # Remove any cells with missing labels..
  complete_cases <- complete.cases(final_labels$cell_type)
  final_labels <- final_labels[complete_cases, ]
  combined_features <- combined_features[final_labels$cell_id, ]
  
  # ----------------------------------------
  # CREATE TRAIN/TEST SPLITS
  # ----------------------------------------
  
  set.seed(42)  # For reproducibility
  
  # Stratified sampling to ensure balanced representation of cell types
  train_indices <- c()
  for (cell_type in unique(final_labels$cell_type)) {
    type_indices <- which(final_labels$cell_type == cell_type)
    n_train <- floor(length(type_indices) * 0.7)  # 70% for training
    train_indices <- c(train_indices, sample(type_indices, n_train))
  }
  
  test_indices <- setdiff(1:nrow(final_labels), train_indices)
  
  # Create split labels
  data_splits <- data.frame(
    cell_id = final_labels$cell_id,
    split = ifelse(1:nrow(final_labels) %in% train_indices, "train", "test")
  )
  
  # ----------------------------------------
  # Export FILES for python :...
  # ----------------------------------------
  
  cat("Exporting Python-ready files...\n")
  
  # Main feature matrices
  write.csv(combined_features, file.path(output_dir_python, "combined_features.csv"))
  write.csv(rna_features_matched, file.path(output_dir_python, "rna_features_2000.csv"))  
  write.csv(atac_features_matched, file.path(output_dir_python, "atac_features_5000.csv"))
  
  # Labels and splits
  write.csv(final_labels, file.path(output_dir_python, "cell_labels.csv"), row.names = FALSE)
  write.csv(data_splits, file.path(output_dir_python, "data_splits.csv"), row.names = FALSE)
  
  # Feature information
  feature_info <- data.frame(
    feature_name = colnames(combined_features),
    feature_type = c(rep("RNA", ncol(rna_features_matched)), 
                     rep("ATAC", ncol(atac_features_matched))),
    feature_index = 1:ncol(combined_features)
  )
  write.csv(feature_info, file.path(output_dir_python, "feature_info.csv"), row.names = FALSE)
  
  # Summary statistics
  summary_stats <- data.frame(
    metric = c("total_cells", "total_features", "rna_features", "atac_features", 
               "train_cells", "test_cells", "cell_types"),
    value = c(nrow(combined_features), ncol(combined_features), ncol(rna_features_matched),
              ncol(atac_features_matched), sum(data_splits$split == "train"),
              sum(data_splits$split == "test"), length(unique(final_labels$cell_type)))
  )
  write.csv(summary_stats, file.path(output_dir_python, "data_summary.csv"), row.names = FALSE)
  
  # Cell type distribution
  cell_type_dist <- table(final_labels$cell_type)
  write.csv(as.data.frame(cell_type_dist), file.path(output_dir_python, "cell_type_distribution.csv"))
  
  cat("=== PYTHON-READY DATA PREPARATION COMPLETED ===\n")
  cat("Files saved in:", output_dir_python, "\n")
  cat("Ready for Python ML pipeline!\n")
  
  return(list(
    features = combined_features,
    labels = final_labels,
    splits = data_splits
  ))
}

# ================================================================================
# MAIN EXECUTION FUNCTIONS....
# ================================================================================

# Function to run Team 1 work
run_team1 <- function() {
  cat("Starting Team 1 RNA processing...\n")
  result <- rna_processing_team1()
  cat("Team 1 completed successfully!\n")
  return(result)
}

# Function to run Team 2 work  
run_team2 <- function() {
  cat("Starting Team 2 ATAC processing...\n")
  result <- atac_processing_team2()
  cat("Team 2 completed successfully!\n")
  return(result)
}

# Function to combine results
run_integration <- function() {
  cat("Starting data integration for Python...\n")
  result <- combine_for_python()
  cat("Integration completed successfully!\n")
  return(result)
}
