#!/usr/bin/env Rscript
# =============================================================================
# Install all R package dependencies for the publication code
# Run once before executing any analysis script.
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Bioconductor packages
bioc_pkgs <- c(
  "csaw", "edgeR", "ChIPseeker", "DESeq2",
  "GenomicAlignments", "GenomicRanges", "GenomicFeatures",
  "Rsamtools", "rtracklayer",
  "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
  "BSgenome.Dmelanogaster.UCSC.dm6", "Biostrings",
  "EBImage", "ComplexHeatmap"
)

# CRAN packages
cran_pkgs <- c(
  "ggplot2", "ggrepel", "pheatmap", "patchwork",
  "dplyr", "tidyr", "tibble", "tidyverse",
  "UpSetR", "ggvenn",
  "tiff", "caret", "tidymodels", "corrr", "ggfortify",
  "GLCMTextures", "reticulate", "umap", "cluster"
)

cat("Installing Bioconductor packages...\n")
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

cat("Installing CRAN packages...\n")
install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

cat("\n✅ All dependencies installed. Verify with sessionInfo().\n")
