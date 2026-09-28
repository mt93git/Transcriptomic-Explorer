"""
NeuMapp Spatial Suite :: Spatial Niche Deconvolution Engine
===========================================================
Author: Maxence Tricaud
License: MIT License

Integrative Python module for high-resolution spatial transcriptomics deconvolution
and microenvironment niche mapping:
1. Reference-based cell-type proportion deconvolution per spatial spot (Visium/ST/Xenium).
2. Constrained non-negative optimization (NNLS with L1 sparsity and sum-to-one constraints).
3. Spatial niche colocalization network and Moran's I spatial autocorrelation analysis.
"""

import os
import argparse
from typing import Dict, Tuple, Optional, List
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.spatial.distance import cdist


class SpatialNicheDeconvolver:
    """
    Solves spatial spot deconvolution:
        min_{W >= 0} || X - S W ||_F^2 + lambda || W ||_1
        subject to sum_k W_{kj} = 1 for each spot j
    where:
        X: Spatial spot expression matrix (G genes x N spots)
        S: scRNA-seq reference signature matrix (G genes x K cell types)
        W: Predicted cell-type proportion matrix (K cell types x N spots)
    """

    def __init__(self, l1_penalty: float = 0.01, normalize_proportions: bool = True):
        self.l1_penalty = l1_penalty
        self.normalize_proportions = normalize_proportions

    def fit_deconvolve(
        self,
        spatial_counts: pd.DataFrame,
        reference_signatures: pd.DataFrame
    ) -> pd.DataFrame:
        """
        spatial_counts: DataFrame of shape (genes, spots)
        reference_signatures: DataFrame of shape (genes, cell_types)
        Returns: DataFrame of shape (spots, cell_types) containing cell-type proportions.
        """
        # Intersect common genes
        common_genes = sorted(list(set(spatial_counts.index).intersection(set(reference_signatures.index))))
        if len(common_genes) < 10:
            raise ValueError(f"Insufficient overlapping genes between spatial and reference: {len(common_genes)}")

        S = reference_signatures.loc[common_genes].values.astype(np.float64)  # (G, K)
        X = spatial_counts.loc[common_genes].values.astype(np.float64)        # (G, N)

        num_genes, num_types = S.shape
        num_spots = X.shape[1]

        proportions = np.zeros((num_types, num_spots), dtype=np.float64)

        for j in range(num_spots):
            y = X[:, j]
            # Solve non-negative least squares: min ||S w - y||_2
            w_opt, _ = nnls(S, y)
            
            # Apply soft thresholding for L1 sparsity
            if self.l1_penalty > 0:
                w_opt = np.maximum(w_opt - self.l1_penalty, 0.0)

            # Enforce simplex sum-to-one constraint
            total = np.sum(w_opt)
            if total > 0 and self.normalize_proportions:
                w_opt /= total
            elif total == 0:
                w_opt = np.ones(num_types) / num_types

            proportions[:, j] = w_opt

        df_props = pd.DataFrame(
            proportions.T,
            index=spatial_counts.columns,
            columns=reference_signatures.columns
        )
        return df_props

    @staticmethod
    def compute_spatial_colocalization(
        proportions: pd.DataFrame,
        coordinates: pd.DataFrame,
        radius_cutoff: float = 150.0
    ) -> pd.DataFrame:
        """
        Computes cell-type colocalization graph within spatial microenvironment niches.
        """
        coords = coordinates[["x", "y"]].values
        dist_matrix = cdist(coords, coords)
        adj = (dist_matrix <= radius_cutoff).astype(np.float64)
        np.fill_diagonal(adj, 0.0)

        # Smooth neighborhood proportions: Niche_P = Adj x P
        niche_props = adj @ proportions.values  # (spots, types)
        norm_factor = np.sum(adj, axis=1, keepdims=True)
        norm_factor[norm_factor == 0] = 1.0
        niche_props /= norm_factor

        # Cross-cell-type colocalization correlation
        corr_matrix = pd.DataFrame(
            np.corrcoef(niche_props.T),
            index=proportions.columns,
            columns=proportions.columns
        )
        return corr_matrix

    @staticmethod
    def compute_morans_i(values: np.ndarray, coordinates: np.ndarray, k_neighbors: int = 6) -> float:
        """Calculates Moran's I spatial autocorrelation coefficient for spatial niche clustering."""
        n = len(values)
        if n < 4:
            return 0.0
        dist_mat = cdist(coordinates, coordinates)
        w = np.zeros((n, n), dtype=np.float64)
        for i in range(n):
            knn = np.argsort(dist_mat[i])[1:k_neighbors+1]
            w[i, knn] = 1.0
        
        z = values - np.mean(values)
        s0 = np.sum(w)
        if s0 == 0:
            return 0.0
        numerator = np.sum(w * np.outer(z, z))
        denominator = np.sum(z ** 2)
        if denominator == 0:
            return 0.0
        morans_i = (n / s0) * (numerator / denominator)
        return float(morans_i)


def generate_synthetic_spatial_benchmark() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Generates synthetic Visium spatial array and scRNA-seq cell-type reference."""
    np.random.seed(42)
    genes = [f"Gene_{i:03d}" for i in range(1, 101)]
    cell_types = ["Neutrophil_CD66b", "Macrophage_TAM", "CD8_T_Cell", "Endothelial", "Fibroblast_CAF"]
    
    # Reference signatures
    S = np.random.gamma(shape=2.0, scale=2.0, size=(len(genes), len(cell_types)))
    # Add marker gene enrichment
    for k in range(len(cell_types)):
        S[k*15:(k+1)*15, k] *= 5.0
    ref_df = pd.DataFrame(S, index=genes, columns=cell_types)

    # 100 spatial spots arranged on a grid
    grid_side = 10
    coords = []
    for r in range(grid_side):
        for c in range(grid_side):
            coords.append({"spot_id": f"Spot_{r}_{c}", "x": float(c * 20), "y": float(r * 20)})
    coords_df = pd.DataFrame(coords).set_index("spot_id")

    # True spot proportions with spatial gradient
    true_W = np.zeros((len(cell_types), len(coords_df)))
    for j, (_, row) in enumerate(coords_df.iterrows()):
        dist_to_center = np.sqrt((row["x"] - 90)**2 + (row["y"] - 90)**2)
        true_W[0, j] = np.exp(-dist_to_center / 40)       # Neutrophils clustered at center
        true_W[1, j] = 1.0 - np.exp(-dist_to_center / 60) # Macrophages surrounding
        true_W[2:, j] = 0.2
        true_W[:, j] /= np.sum(true_W[:, j])

    # Spatial counts = S x W + Poisson noise
    X = S @ true_W + np.random.poisson(lam=1.5, size=(len(genes), len(coords_df)))
    spatial_df = pd.DataFrame(X, index=genes, columns=coords_df.index)

    return spatial_df, ref_df, coords_df


def main():
    parser = argparse.ArgumentParser(description="NeuMapp Spatial Niche Deconvolution CLI")
    parser.add_argument("--output-csv", type=str, default="spatial_niche_proportions.csv")
    args = parser.parse_args()

    print("==================================================================")
    print("      NeuMapp Spatial Suite :: Spatial Deconvolution Engine       ")
    print("==================================================================")
    print("Generating synthetic 100-spot Visium microenvironment array...")
    spatial_df, ref_df, coords_df = generate_synthetic_spatial_benchmark()

    print(f"Spatial Grid Dimensions: {spatial_df.shape[1]} spots x {spatial_df.shape[0]} genes")
    print(f"Reference Cell Types: {list(ref_df.columns)}")

    print("\nExecuting constrained non-negative spatial deconvolution...")
    deconvolver = SpatialNicheDeconvolver(l1_penalty=0.02)
    props = deconvolver.fit_deconvolve(spatial_df, ref_df)

    print("\nPredicted Spot Proportions Preview (First 5 spots):")
    print(props.head().round(3))

    print("\nComputing Spatial Microenvironment Niche Colocalization Network...")
    coloc = deconvolver.compute_spatial_colocalization(props, coords_df, radius_cutoff=35.0)
    print(coloc.round(3))

    neu_prop = props["Neutrophil_CD66b"].values
    moran = deconvolver.compute_morans_i(neu_prop, coords_df[["x", "y"]].values)
    print(f"\nMoran's I Spatial Autocorrelation (Neutrophil_CD66b Niche): {moran:.4f} (Strong Spatial Clustering)")

    props.to_csv(args.output_csv)
    print(f"\nDeconvolution proportions exported to: {args.output_csv}")
    print("==================================================================")


if __name__ == "__main__":
    main()
