# R PROCESSING SCRIPTS - TEAM DIVISION -- 2 teams: 

# ================================================================================
# Team Member 1-Rana Hamed : RNA Data Processing & Cell Type Annotation
# Team Member 2-Syrus+Emmanuel: ATAC Data Processing & Feature Engineering
# Final: Combined Output for Python ML Pipeline....to start the ML building...



# ================================================================================





# Rana Hamed-team 1: RNA DATA PROCESSING & CELL TYPE ANNOTATION ...
# ================================================================================
#first install ll needed packages from the "setup_environment.txt" file...--->
#use this link to download the data :
#https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-10-k-1-standard-1-0-0

#---- the needed files to download ---------------:

#1 filtered_feature_bc_matrix.h5
#2 per_barcode_summary_metrics.csv
#3 atac_fragments.tsv.gz
#4 atac_peaks.bed
#5 atac_peak_annotations.tsv
#6 secondary_analysis.tar.gz
#7 atac_fragment.tsv.gz.tbi

# now ready to only run the following steps... better to run step by step for error easy check : ----
# ================================================================================
# ================================================================================


# load Required Libraries
library(Seurat)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(patchwork)

# ----------------------------------------
# STEP 1: RNA DATA EXTRACTION & QC
# ----------------------------------------
# Setting working directory
setwd("~/singlecell/ATAC")

rna_processing_team1 <- function() {
  
  cat(" RNA Data Processing Started ===\n")
  
  # Setup directories
  input_dir <- "~/singlecell/ATAC"
  output_dir_team1 <- "team1_rna_output"
  if (!dir.exists(output_dir_team1)) dir.create(output_dir_team1, recursive = TRUE)
  
  # ----------------------------------------
  # step1.1: LOAD MAIN DATA FILES..
  # ----------------------------------------
  
  cat("Loading main data matrix...\n")
  
  # Load H5 file (the main data)
  h5_file <- file.path(input_dir, "pbmc_unsorted_10k_filtered_feature_bc_matrix.h5")
  
  if (!file.exists(h5_file)) {
    stop("Main data file not found: ", h5_file)
  }
  
  # Read 10X data
  data_10x <- Read10X_h5(h5_file)
  
  # Extract RNA data (Gene Expression)
  if (is.list(data_10x)) {
    rna_counts <- data_10x$`Gene Expression`
    cat("Multiome data detected - extracted RNA portion\n")
  } else {
    rna_counts <- data_10x
    cat("Single modality data detected\n")
  }
  
  cat("RNA data dimensions:", nrow(rna_counts), "genes x", ncol(rna_counts), "cells\n")
  
  # ----------------------------------------
  # step1.2: LOAD QC METRICS...
  # ----------------------------------------
  
  cat("Loading per-cell QC metrics...\n")
  
  # Load barcode metrics 
  qc_file <- file.path(input_dir, "pbmc_unsorted_10k_per_barcode_metrics.csv")
  
  if (file.exists(qc_file)) {
    qc_metrics <- read.csv(qc_file, row.names = 1)
    cat("Loaded QC metrics for", nrow(qc_metrics), "cells\n")
  } else {
    cat("QC metrics file not found - will calculate basic metrics\n")
    qc_metrics <- NULL
  }
  
  # ----------------------------------------
  # step1.3: CREATE SEURAT OBJECT...with seurat pack-
  # ----------------------------------------
  
  cat("Creating Seurat object...\n")
  
  # Create Seurat object
  pbmc_rna <- CreateSeuratObject(
    counts = rna_counts,
    project = "PBMC_RNA",
    min.cells = 3,      # Filter genes present in <3 cells
    min.features = 200  # Filter cells with <200 genes
  )
  
  # Add QC metrics if available
  if (!is.null(qc_metrics)) {
    # Match cell names and add metrics
    common_cells <- intersect(colnames(pbmc_rna), rownames(qc_metrics))
    pbmc_rna <- subset(pbmc_rna, cells = common_cells)
    
    # Add relevant QC metrics
    if ("num_genes" %in% colnames(qc_metrics)) {
      pbmc_rna$qc_num_genes <- qc_metrics[colnames(pbmc_rna), "num_genes"]
    }
    if ("num_umis" %in% colnames(qc_metrics)) {
      pbmc_rna$qc_num_umis <- qc_metrics[colnames(pbmc_rna), "num_umis"]
    }
  }
  
  # Calculate mitochondrial gene percentage
  pbmc_rna[["percent.mt"]] <- PercentageFeatureSet(pbmc_rna, pattern = "^MT-")
  
  # Calculate ribosomal gene percentage  
  pbmc_rna[["percent.ribo"]] <- PercentageFeatureSet(pbmc_rna, pattern = "^RP[SL]")
  
  cat("Seurat object created with", ncol(pbmc_rna), "cells and", nrow(pbmc_rna), "genes\n")
  
  # ----------------------------------------
  # step 1.4: QUALITY CONTROL FILTERING
  # ----------------------------------------
  
  cat("Applying quality control filters...\n")
  
  # Store pre-filter counts..
  
  n_cells_before <- ncol(pbmc_rna)
  n_genes_before <- nrow(pbmc_rna)
  
  # Define QC thresholds
  min_features <- 500    # Minimum genes per cell...
  max_features <- 7000   # Maximum genes per cell (remove potential doublets)
  min_counts <- 1000     # Minimum UMI counts per cell
  max_mt_percent <- 20   # Maximum mitochondrial percentage
  
  # Apply filters
  pbmc_rna <- subset(pbmc_rna, subset = nFeature_RNA > min_features & 
                       nFeature_RNA < max_features &
                       nCount_RNA > min_counts &
                       percent.mt < max_mt_percent)
  
  # Remove genes expressed in <5 cells after filtering
  counts <- GetAssayData(pbmc_rna, assay = "RNA", slot = "counts")
  genes_to_keep <- rowSums(counts > 0) >= 5
  pbmc_rna <- subset(pbmc_rna, features = rownames(pbmc_rna)[genes_to_keep])
  
  cat("After QC filtering:\n")
  cat("  - Cells:", n_cells_before, "->", ncol(pbmc_rna), 
      "(", round(ncol(pbmc_rna)/n_cells_before*100, 1), "%)\n")
  cat("  - Genes:", n_genes_before, "->", nrow(pbmc_rna), 
      "(", round(nrow(pbmc_rna)/n_genes_before*100, 1), "%)\n")
  
  # ----------------------------------------
  # step1.5: DATA NORMALIZATION & FEATURE SELECTION ....
  # ----------------------------------------
  
  cat("Normalizing data and selecting variable features.....\n")
  
  # Normalize data
  pbmc_rna <- NormalizeData(pbmc_rna, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
  
  # Find variable features
  pbmc_rna <- FindVariableFeatures(pbmc_rna, selection.method = "vst", 
                                   nfeatures = 2000)
  
  # Scale data for PCA
  all_genes <- rownames(pbmc_rna)
  pbmc_rna <- ScaleData(pbmc_rna, features = all_genes)
  
  # ----------------------------------------
  # step1.6: DIMENSIONALITY REDUCTION & CLUSTERING...with PCA and UMAP..
  # ----------------------------------------
  
  cat("Start dimensionality reduction and clustering...\n")
  
  # PCA
  pbmc_rna <- RunPCA(pbmc_rna, features = VariableFeatures(object = pbmc_rna))
  
  # Find neighbors and clusters
  pbmc_rna <- FindNeighbors(pbmc_rna, dims = 1:15)
  pbmc_rna <- FindClusters(pbmc_rna, resolution = 0.5)
  
  #TSNE for visualization
  pbmc_rna <- RunTSNE(pbmc_rna, dims = 1:15)
  
  # UMAP for visualization
  pbmc_rna <- RunUMAP(pbmc_rna, dims = 1:15)
  
  # ----------------------------------------
  # step1.7: CELL TYPE ANNOTATION
  # ----------------------------------------
  
  cat("Performing cell type annotation...\n")
  
  # Find marker genes for each cluster
  rna_markers <- FindAllMarkers(pbmc_rna, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25)
  
  # Automated annotation
  ref = celldex::HumanPrimaryCellAtlasData()
  head(ref)
  
  # Main annotations
  my.sce.new = as.SingleCellExperiment(pbmc_rna)
  preds = SingleR(my.sce.new, ref=ref, labels=ref$label.main)
  table(preds$labels)
  
  # Add cell type annotations to Seurat object 
  pbmc_rna$cell_type <- preds$labels[match(rownames(pbmc_rna@meta.data),rownames(preds))]
  
   cat("Cell type annotation completed:\n")
  print(table(pbmc_rna$cell_type))
  
  # ----------------------------------------
  # step 1.8:  EXPORT RESULTS...:
  # -------------------------------...:---------
  
  cat("Exporting results...\n")
  
  # Extract normalized RNA data for top 2000 variable genes
  variable_genes <- head(VariableFeatures(pbmc_rna), 2000)
  rna_data_normalized <- GetAssayData(pbmc_rna, assay = "RNA", slot = "data")
  rna_features <- as.data.frame(as.matrix(t(rna_data_normalized[variable_genes, ])))
  
  # Export files
  write.csv(rna_features, file.path(output_dir_team1, "rna_features_2000.csv"))
  write.csv(pbmc_rna@meta.data, file.path(output_dir_team1, "rna_metadata.csv"))
  
  # Cell type labels for ML....
  
  cell_labels <- data.frame(
    cell_id = colnames(pbmc_rna),
    cell_type = pbmc_rna$cell_type,
    cluster = pbmc_rna$seurat_clusters
  )
  write.csv(cell_labels, file.path(output_dir_team1, "cell_type_labels.csv"), row.names = FALSE)
  
  # Variable genes list...
  
  write.csv(data.frame(gene = variable_genes), 
            file.path(output_dir_team1, "variable_genes_list.csv"), row.names = FALSE)
  
  # QC plots
  p1 <- VlnPlot(pbmc_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3)
  p2 <- DimPlot(pbmc_rna, reduction = "umap", label = TRUE, group.by = "cell_type")
  
  ggsave(file.path(output_dir_team1, "rna_qc_plots.png"), p1, width = 12, height = 6)
  ggsave(file.path(output_dir_team1, "rna_cell_types_umap.png"), p2, width = 10, height = 8)
  
  # Save Seurat object
  saveRDS(pbmc_rna, file.path(output_dir_team1, "pbmc_rna_processed.rds"))
  
  cat("=== TEAM 1-Syrus RNA PROCESSING COMPLETED ===\n")
  cat("Output files saved in:", output_dir_team1, "\n")
  
  return(pbmc_rna)
}

