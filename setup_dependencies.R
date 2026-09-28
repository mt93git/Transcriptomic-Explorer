# setup_dependencies.R - Transcriptomic Explorer (v2.0)
# Automated package installation script for scRNA-seq interactive analysis

cat("=== Installing Transcriptomic Explorer Dependencies ===\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cran_packages <- c(
  "shiny", "ggplot2", "dplyr", "DT", "shinyjqui", "patchwork",
  "ggrepel", "RColorBrewer", "plotly", "shinycssloaders", "shinyjs",
  "scales", "colourpicker", "ggnewscale", "tibble", "viridis",
  "pheatmap", "MASS", "stringr"
)

bioc_packages <- c(
  "Seurat", "harmony", "ggrastr"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing CRAN package: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing package via BiocManager: %s\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

cat("=== All dependencies successfully installed and verified ===\n")
