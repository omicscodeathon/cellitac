# scATACtf R Dependencies
# R version: 4.5.1

# Installation script for required R packages
# Run this script in R/RStudio to install all necessary packages

# Check R version
cat("Checking R version...\n")
if (getRversion() < "4.5.1") {
  warning("R version 4.5.1 or higher is recommended")
}

# Install Bioconductor if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Core packages for single-cell analysis
cat("\nInstalling Seurat...\n")
install.packages("Seurat")

cat("\nInstalling Signac...\n")
install.packages("Signac")

# Verify installations and versions
cat("\n=== Installed Package Versions ===\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Signac version:", as.character(packageVersion("Signac")), "\n")

# Expected versions:
# Seurat: 5.3.0
# Signac: 1.15.0

cat("\nNote: Required versions are Seurat 5.3.0 and Signac 1.15.0\n")
cat("If versions differ, please update using:\n")
cat("install.packages('Seurat')\n")
cat("install.packages('Signac')\n")
