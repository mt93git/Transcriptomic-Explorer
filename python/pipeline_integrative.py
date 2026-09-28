"""
NeuMapp Explorer :: Integrative Neutrophil Architecture Pipeline
=================================================================
Author: Maxence Tricaud
License: MIT License

Integrates scRNA-seq trajectory inference and spatial tissue hub deconvolution
for neutrophil developmental ontogeny and tissue compartmentalization:
1. Reconstructs differentiation trajectories, chronomic aging, and pseudotime vectors from single-cell profiles.
2. Derives dynamic signature weights for evolving neutrophil maturation states.
3. Deconvolves 2D spatial microenvironment arrays and tissue compartments to project
   neutrophil subsets into physical anatomical hubs (Bone Marrow, Blood, Marginal Pools, Target Organs).
"""

import sys
import os
import numpy as np
import pandas as pd

from .trajectory_inference import SingleCellTrajectoryEngine, generate_synthetic_trajectory_data
from .spatial_deconvolution import SpatialNicheDeconvolver, generate_synthetic_spatial_benchmark


def run_integrative_pipeline():
    print("==========================================================================")
    print("  NeuMapp Explorer :: Neutrophil Architecture & Spatial Hub Pipeline     ")
    print("==========================================================================")

    # 1. Single-Cell Trajectory Inference
    print("\n--- Step 1: Inferring Maturation & Chronomics Trajectory from scRNA-seq ---")
    df_expr, cluster_labels, _ = generate_synthetic_trajectory_data(n_cells=120)
    traj_engine = SingleCellTrajectoryEngine(n_neighbors=12, n_components=3)
    diff_coords, evals, T_markov = traj_engine.compute_diffusion_map(df_expr.values)
    pseudotime = traj_engine.compute_pseudotime(T_markov, root_cell_idx=0)
    print(f"Successfully computed Diffusion Pseudotime for {len(pseudotime)} cells.")
    print(f"Diffusion Components extracted: {diff_coords.shape[1]}")

    # 2. Spatial Niche Deconvolution
    print("\n--- Step 2: Deconvolving Microenvironment Niches from Spatial Data ---")
    spatial_df, ref_df, coords_df = generate_synthetic_spatial_benchmark()
    deconvolver = SpatialNicheDeconvolver(l1_penalty=0.01)
    props = deconvolver.fit_deconvolve(spatial_df, ref_df)
    print(f"Successfully deconvolved {props.shape[0]} spatial spots across {props.shape[1]} cell types.")

    # 3. Spatial Niche Colocalization
    print("\n--- Step 3: Computing Spatial Colocalization Graph ---")
    coloc = deconvolver.compute_spatial_colocalization(props, coords_df, radius_cutoff=35.0)
    print("Spatial Niche Cross-Correlation Matrix:")
    print(coloc.round(3))

    print("\n==========================================================================")
    print(">>> INTEGRATIVE MULTI-MODAL ANALYSIS SUCCESSFULLY COMPLETED! <<<")
    print("==========================================================================")


if __name__ == "__main__":
    run_integrative_pipeline()
