# NeuMapp Explorer (Transcriptomic-Explorer)

[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![Python Version](https://img.shields.io/badge/Python-%3E%3D3.9-blue.svg)](https://www.python.org/)
[![Shiny](https://img.shields.io/badge/Framework-Shiny%20%7C%20bslib-2C3E50.svg)](https://shiny.posit.co/)
[![Bioconductor](https://img.shields.io/badge/Bioc-Seurat%20v5-brightgreen.svg)](https://satijalab.org/seurat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Interactive Multi-Omic Platform for Exploring the Transcriptomic Architecture, Trajectory Dynamics, and Spatial Tissue Compartmentalization of Neutrophil Hubs.**  
> **Lead Architect & Developer:** Maxence Tricaud (`mtricaud.cetri@gmail.com`)

---

## 1. Executive Summary & Scientific Scope

**NeuMapp Explorer** (maintained in this repository under `Transcriptomic-Explorer`) is an interactive multi-omic computational platform engineered for investigating the transcriptomic landscape, developmental trajectory dynamics, and tissue compartmentalization architecture of neutrophil functional states. The suite bridges single-cell transcriptomic resolution (scRNA-seq) with tissue-level anatomical hubs, combining an interactive **R/Shiny** visualization frontend with an optimized scientific **Python** computing backend for trajectory inference and spatial niche deconvolution.

```
                      ┌─────────────────────────────────────────────────────────────┐
                      │              NeuMapp Explorer Architecture                  │
                      └──────────────────────────────┬──────────────────────────────┘
                                                     │
                     ┌───────────────────────────────┴───────────────────────────────┐
                     ▼                                                               ▼
        ┌─────────────────────────┐                                     ┌─────────────────────────┐
        │  R/Shiny User Frontend  │                                     │  Python Analytical Core │
        ├─────────────────────────┤                                     ├─────────────────────────┤
        │ • modules/neumapp_      │                                     │ • spatial_              │
        │   spatial/ (3400+ loc)  │                                     │   deconvolution.py      │
        │ • modules/seurat_       │◄────────── Cross-Modality ──────────►│ • trajectory_           │
        │   singlecell/           │             Coordination            │   inference.py          │
        │ • launch_suite.R        │                                     │ • pipeline_             │
        │ • root app.R (Seurat v5)│                                     │   integrative.py        │
        └─────────────────────────┘                                     └─────────────────────────┘
```

---

## 2. Biological Paradigm: Neutrophil Spatial Architecture & Compartmentalization

### Resolving the Concept of "Spatial" in Neutrophil Biology
In neutrophil immunobiology, **"spatial" architecture does not merely denote 2D pixel coordinates on an artificial glass slide**, but fundamentally refers to the **anatomical compartmentalization and organ-specific vascular hubs** pioneered by **Andrés Hidalgo et al.** (*Cell* 2019, *Cell* 2020, *Nature Immunology* 2021):

1. **Bone Marrow Precursor & Retention Niche:** Pro-Neu and pre-Neu proliferative reservoirs (`c-Kit+`, `CXCR4+`).
2. **Circulating Vascular Pool:** Diurnally oscillating, mature effector neutrophils (`CD62L_high`, `CXCR2+`).
3. **Marginal Vascular Pools:** Non-canonical intravascular marginated reservoirs adhering in lung capillary beds and splenic red pulp.
4. **Target Tissue & Microenvironment Hubs:** Extravasated neutrophil phenotypes adapting to localized niches, including acute inflammatory foci and tumor microenvironments (TANs: N1 anti-tumor vs N2 immunosuppressive niches).

### The Computational Challenge Addressed by NeuMapp Explorer
* **Continuous Lineage Kinetics (scRNA-seq):** Neutrophil maturation and diurnal chronomic aging represent a continuous phenotypic continuum rather than discrete cell states. NeuMapp resolves this using **Spectral Diffusion Maps** and **Diffusion Pseudotime (DPT)**.
* **Tissue Compartmentalization Deconvolution:** Bulk tissue biopsies, spatial arrays (10x Visium/Xenium), and organ marginal pools represent mixtures of multiple cellular niches. NeuMapp resolves the exact fractional distribution of neutrophil functional states via **Constrained Non-Negative Matrix Factorization (NNLS + $\ell_1$ sparsity)** and spatial autocorrelation (**Moran's $I$**).

---

## 3. Algorithmic Modules & Mathematical Foundations

### 3.1 Spatial Niche & Tissue Compartmentalization Engine (`python/spatial_deconvolution.py`)
Deconvolves mixed spot or tissue compartment profiles $\mathbf{X} \in \mathbb{R}^{G \times N}$ against single-cell reference signature matrices $\mathbf{S} \in \mathbb{R}^{G \times K}$:

$$\min_{\mathbf{W} \ge 0} \|\mathbf{X} - \mathbf{S}\mathbf{W}\|_F^2 + \lambda \|\mathbf{W}\|_1 \quad \text{subject to} \quad \sum_{k=1}^K W_{kj} = 1 \quad \forall j$$

* **Spatial Autocorrelation (Moran's $I$):** Identifies non-random spatial clustering of specific neutrophil hubs across tissue coordinates:
  $$I = \frac{N}{S_0} \frac{\sum_i \sum_j w_{ij}(z_i - \bar{z})(z_j - \bar{z})}{\sum_i (z_i - \bar{z})^2}$$
* **Niche Colocalization Graph:** Generates a cross-correlation network quantifying spatial co-occurrence or mutual exclusion between distinct neutrophil states and surrounding stromal/immune cell types.

### 3.2 Single-Cell Trajectory & Chronomics Engine (`python/trajectory_inference.py`)
Reconstructs non-linear developmental and circadian aging paths along the single-cell manifold:
* **Adaptive Gaussian Transition Kernel:**
  $$K(x_i, x_j) = \exp\left(-\frac{\|x_i - x_j\|^2}{\sigma(x_i)\sigma(x_j)}\right)$$
  where $\sigma(x_i)$ is determined by the distance to the $k$-th nearest neighbor.
* **Spectral Diffusion Map:** Normalizes kernel into a row-stochastic Markov transition operator $\mathbf{T}$, computing top diffusion components (DC1, DC2, DC3) capturing global lineage topology.
* **Diffusion Pseudotime (DPT):** Measures random-walk geodesic distances from the bone marrow precursor root (`Pre_Neu`) to terminal chronomic and activation states.
* **PAGA-Style Markovian Connectivity:** Computes coarse-grained cluster transition probabilities between maturation stages.

### 3.3 Integrative Cross-Modality Runner (`python/pipeline_integrative.py`)
Couples single-cell trajectory kinetics with spatial compartment deconvolution, mapping continuous developmental vectors directly into physical organ microenvironments.

### 3.4 Interactive R/Shiny Dashboards (`modules/`)
* **`modules/neumapp_spatial/` (3,400+ lines):** Histology slice overlay, spatial feature contours, spot clustering, and differential microenvironment niche exploration.
* **`modules/seurat_singlecell/` & root `app.R` (1,360+ lines):** Native support for Seurat v5 Assay5 multi-layer architectures, real-time dynamic subsetting, on-the-fly UMAP re-clustering, cell density glow filters, and automated differential marker detection.

---

## 4. ⚡ Quick Start & Verification

### Option A: Python Analytical Core (Trajectory & Deconvolution)

```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# 2. Run the end-to-end integrative pipeline (direct execution or module mode)
python3 python/pipeline_integrative.py
# or: python3 -m python.pipeline_integrative

# 3. Run standalone neutrophil spatial hub deconvolution CLI
python3 python/spatial_deconvolution.py --output-csv spatial_niche_proportions.csv

# 4. Run standalone single-cell trajectory & chronomics CLI
python3 python/trajectory_inference.py --n-neighbors 15 --output-csv trajectory_pseudotime_results.csv
```

---

### Option B: Interactive R/Shiny Dashboards

```bash
# 1. Install R dependencies
Rscript setup_dependencies.R

# 2. Launch root Seurat v5 single-cell explorer
Rscript -e "shiny::runApp('app.R')"
```
*Load the included demo dataset (`demo_data/demo_pbmc_small.rds`) directly in the UI for instant testing.*

To launch the multi-module suite selector:
```R
source("launch_suite.R")
launch_suite("spatial")      # Launches NeuMapp Spatial Suite
# OR
launch_suite("singlecell")   # Launches Seurat Single-Cell Explorer
```

---

## 5. Repository Structure

```
Transcriptomic-Explorer/
├── launch_suite.R              # Interactive multi-module launcher
├── app.R                       # Single-cell Shiny dashboard (root, Seurat v5)
├── setup_dependencies.R        # Automated CRAN/Bioconductor dependency installer
├── pyproject.toml              # PEP 621 Python packaging configuration
├── requirements.txt            # Python dependencies (numpy, scipy, pandas)
├── python/                     # NeuMapp Python Analytical Engine
│   ├── __init__.py             # Package exports
│   ├── spatial_deconvolution.py# Constrained NNLS & Moran's I spatial niche engine
│   ├── trajectory_inference.py # Spectral diffusion maps & Markov transition graph
│   └── pipeline_integrative.py # Integrative scRNA-seq trajectory + spatial hub runner
├── modules/
│   ├── neumapp_spatial/        # Spatial Transcriptomics Shiny Suite (3400+ lines)
│   │   ├── server.R
│   │   └── ui.R
│   └── seurat_singlecell/      # Seurat v5 Manifold Explorer
│       └── app.R
├── data/
│   └── sample_metadata.csv     # Sample and cohort metadata schema
├── demo_data/
│   └── demo_pbmc_small.rds     # Curated demo dataset for instant evaluation
├── LICENSE                     # MIT Open-Source License
├── CITATION.cff                # Academic citation metadata
├── CONTRIBUTING.md             # Developer guidelines
├── SECURITY.md                 # Security and vulnerability reporting
└── README.md
```

---

## 6. Scientific References & Conceptual Foundation

The architectural design of **NeuMapp Explorer** is founded on the following landmark works on neutrophil compartmentalization, vascular hubs, and heterogeneity:

1. **Hidalgo, A., et al.** (2019). *Neutrophils: Forging the Future of Immunology.* **Cell**, 179(3), 585–597.
2. **Ballesteros, I., et al.** (2020). *Cellular and Molecular Determinants of Neutrophil Heterogeneity and Aging across Tissues.* **Cell**, 181(4), 842–859.
3. **Casanova-Acebes, M., et al.** (2018). *Neutrophils Instruct Homeostatic and Pathological States in Distinct Organ Niches.* **Cell**, 174(5), 1170–1182.
4. **Hidalgo, A., et al.** (2021). *Functional Epigenomics and Vascular Niches of Leukocyte Margination.* **Nature Immunology**, 22(8), 940–953.

---

## 7. License & Citation

Distributed under the **MIT License**. Copyright © 2025–2026 Maxence Tricaud.

```bibtex
@software{tricaud2026neumapp,
  author       = {Tricaud, Maxence},
  title        = {NeuMapp Explorer: Multi-Omic Landscape of Neutrophil Architecture, Trajectory Dynamics, and Spatial Tissue Compartmentalization},
  year         = {2026},
  url          = {https://github.com/mt93git/Transcriptomic-Explorer},
  version      = {2.1.0}
}
```
