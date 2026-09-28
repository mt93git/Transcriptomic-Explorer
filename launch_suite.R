#!/usr/bin/env Rscript
# =============================================================================
# NeuMapp Spatial Suite :: Unified Multi-Omics Launcher
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
module <- if (length(args) > 0) args[1] else "1"

cat("====================================================================\n")
cat("            NeuMapp Spatial Suite :: Application Launcher           \n")
cat("====================================================================\n")
cat("Available Platforms:\n")
cat("  1. NeuMapp Spatial Explorer (modules/neumapp_spatial/)\n")
cat("  2. Seurat scRNA-seq Manifold Explorer (modules/seurat_singlecell/)\n")
cat("====================================================================\n\n")

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The 'shiny' package is required. Install with: install.packages('shiny')")
}

if (module == "1") {
  cat("Launching NeuMapp Spatial Explorer ...\n")
  shiny::runApp("modules/neumapp_spatial")
} else if (module == "2") {
  cat("Launching Seurat Single-Cell Manifold Explorer ...\n")
  shiny::runApp("modules/seurat_singlecell/app.R")
} else {
  cat("Launching default NeuMapp Spatial Explorer ...\n")
  shiny::runApp("modules/neumapp_spatial")
}
