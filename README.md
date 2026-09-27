# Transcriptomic Explorer (v2.0) — NeuMapp Initiative

[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Framework-Shiny%20%7C%20bslib-2C3E50.svg)](https://shiny.posit.co/)
[![Bioconductor](https://img.shields.io/badge/Bioc-Seurat%20v5-brightgreen.svg)](https://satijalab.org/seurat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Interactive Single-Cell RNA-seq Analytical Dashboard & Exploratory Frontend of the NeuMapp Suite.**  
> *Developed by Maxence Tricaud (Libreros Lab — Yale University / Université Laval).*

---

## 1. Overview & Scope within NeuMapp

**Transcriptomic Explorer** is the dedicated interactive exploration engine of the **NeuMapp** computational initiative. Built on top of the **Seurat** ecosystem (with native support for Seurat v5 Assay5 multi-layer architectures), it allows biomedical researchers to interrogate cell state heterogeneity, inspect high-resolution manifold embeddings, and compute differential gene expression without requiring manual coding.

---

## 2. Key Analytical Features

* **Advanced Manifold Visualization:** Real-time 2D UMAP projection with density contours, cell density glow filters, and multi-condition split facet modes.
* **Dynamic On-the-Fly Subsetting & Reprojection:** Graphically select cell clusters or metadata subsets and recalculate normalization, scaling, PCA, and UMAP embeddings in real time.
* **Marker Discovery Backend:** Automated differential expression testing (One-vs-All or Pairwise contrasts) coupled to hierarchical heatmap generation.
* **Vector Graphics Export:** Direct download of publication-ready PDF and SVG plots.

---

## 3. ⚡ Quick Start & Verification

### Launch in RStudio
```R
# 1. Clone the repository
# git clone https://github.com/mt93git/Transcriptomic-Explorer.git

# 2. Launch the application
shiny::runApp("app.R")
```

### Verification with Built-in Seurat Data
To immediately test the interface without an external `.rds` file, run this snippet in R:
```R
library(Seurat)
# Save Seurat's built-in PBMC dataset for immediate upload
saveRDS(pbmc_small, "demo_pbmc.rds")
# Launch app and select demo_pbmc.rds in the file uploader
shiny::runApp("app.R")
```

---

## 4. Repository Structure

```
Transcriptomic-Explorer/
├── app.R                       # Full interactive Shiny dashboard (1300+ lines)
├── LICENSE                     # MIT Open-Source License
├── .gitignore
└── README.md
```

---

## 5. License & Citation

Distributed under the **MIT License**. See `LICENSE` for details.

Citation:
> Tricaud M. *Transcriptomic Explorer: Interactive Single-Cell RNA-seq Exploration Engine for the NeuMapp Initiative*. Yale University / Université Laval (2025-2026). GitHub: `https://github.com/mt93git/Transcriptomic-Explorer`.
