# Transcriptomic Explorer (v2.0)

**An R/Shiny Dashboard for Interactive Single-Cell Analysis**

## 1. Overview
The **Transcriptomic Explorer** is a lightweight, standalone application designed to facilitate the exploration of Single-Cell RNA-seq (scRNA-seq) data. It provides a graphical interface for the **Seurat** computational framework, allowing researchers to visualize dimensionality reduction, perform real-time subsetting, and calculate differential expression without writing code.

**Primary Capabilities:**
* **Visualization:** Interactive UMAP plotting with customizable density, contour, and "glow" overlays.
* **Comparative Analysis:** "Split-View" and "Side-by-Side" modes for comparing experimental conditions.
* **Subsetting & Reprojection:** Select cell populations based on metadata or gene expression thresholds and re-run the analysis pipeline (Normalize ? Scale ? PCA ? UMAP) in real-time.
* **Differential Expression:** Automated identification of marker genes (One-vs-All or Pairwise comparisons) with integrated Heatmaps.

## 2. Key Features

### ?? Dimensionality Reduction (UMAP)
* **Sided/Split Visualization:** Facet plots by metadata (e.g., \Condition\, \Sample\) to assess batch effects or biological shifts.
* **Density Overlays:** Toggle density contours and greyscale modes to visualize cell concentration in high-density clusters.
* **Customizable Aesthetics:** Granular control over point size, transparency, and color palettes (Viridis, RColorBrewer, Custom Hex).

### ?? Real-Time Analysis Pipeline
* **Dynamic Re-Clustering:** Adjust Leiden/Louvain resolution sliders to explore clustering granularity instantly.
* **Advanced Subsetting Engine:** Filter cells using:
    * **Categorical Metadata:** (e.g., specific clusters or cell types).
    * **Numerical Features:** Filter based on specific gene expression levels.
    * **Composed Markers:** Compute "Signature Scores" (e.g., Inflammation Score, Cell Cycle Score) and filter cells that exceed a specific threshold.
* **Reprojection Logic:** When subsetting, the app intelligently detects **Seurat v5** vs **Legacy** objects. It re-processes the subset (JoinLayers if necessary, Normalize, Scale, PCA, UMAP) to generate an accurate projection of the isolated population.

### ?? Marker & Distribution Analysis
* **Multi-Gene Visualization:** Dot Plots and Feature Grids for visualizing gene lists across groups.
* **Composition Analysis:** Stacked bar charts and donut plots to visualize cluster composition across experimental groups.
* **Differential Expression (DEG):**
    * **Pairwise:** Compare Group A vs. Group B.
    * **Global:** Identify markers distinguishing one cluster from all others.
    * **Export:** Download DEG tables as CSV.

## 3. Requirements & Dependencies

The application requires **R (v4.0+)** and the following packages:

* **Core:** \shiny\, \Seurat\ (v4 or v5 compatible), \dplyr\
* **Visualization:** \ggplot2\, \plotly\, \pheatmap\, \ggrepel\, \RColorBrewer\, \iridis\, \colourpicker\, \ggrastr\
* **UI/UX:** \shinyjqui\ (resizable plots), \shinycssloaders\, \shinyjs\, \DT\

## 4. Usage Instructions

### Step 1: Setup
Ensure all dependencies are installed. You can run the included initialization block in the script, or install manually:
\\\
install.packages(c("shiny", "Seurat", "ggplot2", "dplyr", "DT", "shinyjqui", "patchwork", "ggrepel", "RColorBrewer", "plotly", "shinycssloaders", "shinyjs", "scales", "colourpicker", "ggnewscale", "tibble", "viridis", "pheatmap", "harmony", "ggrastr", "MASS", "stringr"))
\\\

### Step 2: Launch
Run the application from RStudio:
\\\
shiny::runApp("app.R")
\\\

### Step 3: Load Data
* Upload a pre-processed \.rds\ file (Seurat Object).
* **Note:** The object must contain a **UMAP** dimensionality reduction slot.
* **File Size Limit:** Default is set to 3GB.

## 5. Version History
* **v2.0:** Added Sided/Split UMAP views, Density overlays, Composed Marker Scoring, and automated Seurat v5 compatibility patches (Assay5 JoinLayers).
* **v1.0:** Initial release with basic visualization and DEG features.

## 6. Citation
If used in research, please cite the Seurat R package (Satija Lab) and this repository.
