# cellitac R Dependencies
# R version: 4.4.3 or higher
# Installation script for required R packages
# Run this script in R/RStudio to install all necessary packages

# Check R version
cat("Checking R version...\n")
if (getRversion() < "4.4.3") {
  warning("R version 4.4.3 or higher is required")
}

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# CRAN packages
cat("\nInstalling Seurat...\n")
install.packages("Seurat")

cat("\nInstalling Signac...\n")
install.packages("Signac")

cat("\nInstalling hdf5r...\n")
install.packages("hdf5r")

cat("\nInstalling Matrix...\n")
install.packages("Matrix")

# Bioconductor packages
cat("\nInstalling SingleR...\n")
BiocManager::install("SingleR")

cat("\nInstalling celldex...\n")
BiocManager::install("celldex")

cat("\nInstalling EnsDb.Hsapiens.v75...\n")
BiocManager::install("EnsDb.Hsapiens.v75")

# Verify installations
cat("\n=== Installed Package Versions ===\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Signac version:", as.character(packageVersion("Signac")), "\n")
cat("hdf5r version:", as.character(packageVersion("hdf5r")), "\n")
cat("Matrix version:", as.character(packageVersion("Matrix")), "\n")
cat("SingleR version:", as.character(packageVersion("SingleR")), "\n")
cat("celldex version:", as.character(packageVersion("celldex")), "\n")
cat("EnsDb.Hsapiens.v75 version:", as.character(packageVersion("EnsDb.Hsapiens.v75")), "\n")

cat("\nAll R packages installed successfully!\n")
