"""
NeuMapp Spatial & Trajectory Analytics Suite (Python Engine)
Author: Maxence Tricaud
License: MIT
"""

from .spatial_deconvolution import (
    SpatialNicheDeconvolver,
    generate_synthetic_spatial_benchmark,
)
from .trajectory_inference import (
    SingleCellTrajectoryEngine,
    generate_synthetic_trajectory_data,
)
from .pipeline_integrative import (
    run_integrative_pipeline,
)

__version__ = "2.1.0"
__author__ = "Maxence Tricaud"
__all__ = [
    "SpatialNicheDeconvolver",
    "generate_synthetic_spatial_benchmark",
    "SingleCellTrajectoryEngine",
    "generate_synthetic_trajectory_data",
    "run_integrative_pipeline",
]
