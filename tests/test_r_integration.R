#!/usr/bin/env Rscript
# =============================================================================
# Unit & Integration Tests: NeuMapp Explorer R Modules & Launchers
# =============================================================================
# Author: Maxence Tricaud
# License: MIT
# =============================================================================

cat("====================================================================\n")
cat("      NeuMapp Explorer :: R Architecture & Integration Test Suite    \n")
cat("====================================================================\n")

# 1. Verify launcher function
cat("\n[1/4] Testing launch_suite.R function export ...\n")
stopifnot(file.exists("launch_suite.R"))
source("launch_suite.R")
stopifnot(exists("launch_suite"), is.function(launch_suite))
formals_names <- names(formals(launch_suite))
stopifnot("mode" %in% formals_names)
cat("  ✔ launch_suite() function successfully declared with valid arguments.\n")

# 2. Verify Shiny Spatial Module Architecture
cat("\n[2/4] Testing modules/neumapp_spatial structure ...\n")
spatial_ui <- "modules/neumapp_spatial/ui.R"
spatial_server <- "modules/neumapp_spatial/server.R"
stopifnot(file.exists(spatial_ui), file.exists(spatial_server))
cat(sprintf("  ✔ Spatial UI verified (%d bytes)\n", file.info(spatial_ui)$size))
cat(sprintf("  ✔ Spatial Server verified (%d bytes)\n", file.info(spatial_server)$size))

# 3. Verify Single-Cell Explorer Module
cat("\n[3/4] Testing modules/seurat_singlecell and root app.R ...\n")
sc_app <- "modules/seurat_singlecell/app.R"
root_app <- "app.R"
stopifnot(file.exists(sc_app), file.exists(root_app))
cat("  ✔ Seurat single-cell apps present at module and root levels.\n")

# 4. Verify Demo Dataset Ingestion
cat("\n[4/4] Testing demo dataset integrity ...\n")
demo_rds <- "demo_data/demo_pbmc_small.rds"
stopifnot(file.exists(demo_rds))
pbmc_data <- readRDS(demo_rds)
stopifnot(!is.null(pbmc_data))
cat(sprintf("  ✔ Demo dataset loaded successfully: class '%s'.\n", class(pbmc_data)[1]))

cat("\n====================================================================\n")
cat(">>> ALL R ARCHITECTURE AND INTEGRATION TESTS PASSED (100%) <<<\n")
cat("====================================================================\n")
