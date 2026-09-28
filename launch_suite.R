#!/usr/bin/env Rscript
# =============================================================================
# NeuMapp Explorer :: Unified Multi-Omics Application Launcher
# =============================================================================

cat("====================================================================\n")
cat("            NeuMapp Explorer :: Application Launcher                \n")
cat("====================================================================\n")
cat("Available Platforms:\n")
cat("  - 'spatial'    : NeuMapp Spatial Explorer (modules/neumapp_spatial/)\n")
cat("  - 'singlecell' : Seurat scRNA-seq Manifold Explorer (modules/seurat_singlecell/)\n")
cat("Usage in R session:\n")
cat("  source('launch_suite.R')\n")
cat("  launch_suite('spatial')   # or launch_suite('singlecell')\n")
cat("====================================================================\n\n")

launch_suite <- function(mode = c("spatial", "singlecell"), port = 3838, launch.browser = TRUE, ...) {
  if (missing(mode) || is.null(mode)) {
    mode <- "spatial"
  }
  
  if (is.numeric(mode)) {
    mode <- if (mode == 2) "singlecell" else "spatial"
  } else {
    mode <- match.arg(tolower(as.character(mode)), c("spatial", "singlecell"))
  }
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The 'shiny' package is required. Install with: install.packages('shiny')")
  }
  
  # Resolve base directory relative to this script
  base_dir <- getwd()
  if (!dir.exists(file.path(base_dir, "modules"))) {
    # Check parent directory
    if (dir.exists(file.path(dirname(base_dir), "modules"))) {
      base_dir <- dirname(base_dir)
    }
  }
  
  if (mode == "spatial") {
    target_dir <- file.path(base_dir, "modules", "neumapp_spatial")
    cat(sprintf("Launching NeuMapp Spatial Explorer from: %s\n", target_dir))
    shiny::runApp(target_dir, port = port, launch.browser = launch.browser, ...)
  } else {
    target_file <- file.path(base_dir, "modules", "seurat_singlecell", "app.R")
    cat(sprintf("Launching Seurat Single-Cell Manifold Explorer from: %s\n", target_file))
    shiny::runApp(target_file, port = port, launch.browser = launch.browser, ...)
  }
}

# Execute automatically if run directly via Rscript in non-interactive terminal
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    arg_mode <- tolower(args[1])
    target <- if (arg_mode %in% c("2", "singlecell", "single_cell", "seurat")) "singlecell" else "spatial"
    launch_suite(target)
  }
}
