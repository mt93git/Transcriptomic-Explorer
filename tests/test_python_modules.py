#!/usr/bin/env python3
"""
Unit and Integration Tests for NeuMapp Explorer Python Engine
============================================================
Author: Maxence Tricaud
License: MIT
"""

import sys
import os
import unittest
import numpy as np
import pandas as pd

# Add repo root and python dir to sys.path
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

python_dir = os.path.join(repo_root, "python")
if python_dir not in sys.path:
    sys.path.insert(0, python_dir)

from python.spatial_deconvolution import (
    SpatialNicheDeconvolver,
    generate_synthetic_spatial_benchmark
)
from python.trajectory_inference import (
    SingleCellTrajectoryEngine,
    generate_synthetic_trajectory_data
)
from python.pipeline_integrative import run_integrative_pipeline


class TestNeuMappSpatialDeconvolution(unittest.TestCase):
    """Tests constrained NNLS deconvolution and spatial autocorrelation."""

    def setUp(self):
        self.spatial_df, self.ref_df, self.coords_df = generate_synthetic_spatial_benchmark()
        self.deconvolver = SpatialNicheDeconvolver(l1_penalty=0.01)

    def test_synthetic_data_dimensions(self):
        self.assertEqual(self.spatial_df.shape[0], 100)  # 100 genes
        self.assertEqual(self.spatial_df.shape[1], 100)  # 100 spots
        self.assertEqual(self.ref_df.shape[1], 5)       # 5 cell types / hubs

    def test_deconvolution_proportions_sum_to_one(self):
        props = self.deconvolver.fit_deconvolve(self.spatial_df, self.ref_df)
        self.assertEqual(props.shape[0], 100)
        self.assertEqual(props.shape[1], 5)
        # Check sum to 1 per spot
        row_sums = props.sum(axis=1).values
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-3)
        # Check non-negativity
        self.assertTrue((props.values >= 0).all())

    def test_morans_i_calculation(self):
        props = self.deconvolver.fit_deconvolve(self.spatial_df, self.ref_df)
        neu_prop = props.iloc[:, 0].values
        moran = self.deconvolver.compute_morans_i(neu_prop, self.coords_df[["x", "y"]].values)
        self.assertIsInstance(moran, float)
        self.assertTrue(-1.0 <= moran <= 1.0)

    def test_spatial_colocalization(self):
        props = self.deconvolver.fit_deconvolve(self.spatial_df, self.ref_df)
        coloc = self.deconvolver.compute_spatial_colocalization(props, self.coords_df, radius_cutoff=35.0)
        self.assertEqual(coloc.shape, (5, 5))
        # Diagonal should be 1.0
        np.testing.assert_allclose(np.diag(coloc.values), 1.0, atol=1e-3)


class TestNeuMappTrajectoryInference(unittest.TestCase):
    """Tests spectral diffusion maps, pseudotime, and PAGA connectivity."""

    def setUp(self):
        self.df_expr, self.labels, self.df_truth = generate_synthetic_trajectory_data(n_cells=80)
        self.traj_engine = SingleCellTrajectoryEngine(n_neighbors=10, n_components=3)

    def test_diffusion_map_computation(self):
        diff_coords, evals, T_markov = self.traj_engine.compute_diffusion_map(self.df_expr.values)
        self.assertEqual(diff_coords.shape, (80, 3))
        self.assertEqual(len(evals), 3)
        self.assertEqual(T_markov.shape, (80, 80))
        # Markov transition operator rows should sum to 1
        np.testing.assert_allclose(T_markov.sum(axis=1), 1.0, atol=1e-3)

    def test_diffusion_pseudotime(self):
        _, _, T_markov = self.traj_engine.compute_diffusion_map(self.df_expr.values)
        pseudotime = self.traj_engine.compute_pseudotime(T_markov, root_cell_idx=0)
        self.assertEqual(len(pseudotime), 80)
        self.assertAlmostEqual(pseudotime.min(), 0.0, places=3)
        self.assertAlmostEqual(pseudotime.max(), 1.0, places=3)

    def test_paga_connectivity(self):
        _, _, T_markov = self.traj_engine.compute_diffusion_map(self.df_expr.values)
        df_paga = self.traj_engine.compute_paga_connectivity(T_markov, self.labels)
        unique_labels = len(np.unique(self.labels))
        self.assertEqual(df_paga.shape, (unique_labels, unique_labels))


class TestIntegrativePipeline(unittest.TestCase):
    """Tests end-to-end integrative runner."""

    def test_run_integrative_pipeline_executes(self):
        # Should execute cleanly without throwing
        try:
            run_integrative_pipeline()
        except Exception as e:
            self.fail(f"run_integrative_pipeline raised unexpected exception: {e}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
