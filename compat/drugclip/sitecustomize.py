"""Narrow compatibility allowlist for the official legacy DrugCLIP checkpoint.

This module is injected only into the DrugCLIP retrieval subprocess. It keeps
PyTorch's weights-only loader enabled and allowlists the exact harmless Python
and NumPy metadata types reported by the pinned official checkpoint.
"""

import argparse

import numpy as np
import torch


torch.serialization.add_safe_globals(
    [
        argparse.Namespace,
        np.core.multiarray.scalar,
        np.dtype,
        type(np.dtype(np.float32)),
        type(np.dtype(np.float64)),
        type(np.dtype(np.int64)),
    ]
)

