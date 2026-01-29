import joblib
import numpy as np
from pathlib import Path
from typing import Tuple, Optional


class UMAPProjector:
    """
    Wrapper around a pre-trained UMAP model.
    Responsible ONLY for projection (transform), never fit.
    """

    def __init__(self, model_path: Optional[str] = None):
        """
        Parameters
        ----------
        model_path : str, optional
            Path to serialized UMAP model (.pkl)
            If None, model must be loaded later via load()
        """
        self.umap = None
        if model_path is not None:
            self.load(model_path)

    def load(self, model_path: str):
        """Load UMAP model from file."""
        path = Path(model_path)
        if not path.exists():
            raise FileNotFoundError(f"UMAP model not found: {model_path}")
        try:
            self.umap = joblib.load(model_path)
        except Exception as e:
            raise RuntimeError(f"Failed to load UMAP model from {model_path}") from e

    def transform(self, features: np.ndarray) -> Tuple[float, float]:
        """
        Project a single molecule feature vector into UMAP space.

        Parameters
        ----------
        features : np.ndarray
            Shape (1, n_features) or (n_features,)

        Returns
        -------
        (x, y) : Tuple[float, float]
            2D UMAP coordinates
        """
        if self.umap is None:
            raise RuntimeError("UMAP model not loaded. Call load() first.")

        # Ensure 2D input
        if features.ndim == 1:
            features = features.reshape(1, -1)
        elif features.ndim == 2 and features.shape[0] == 1:
            pass  # Already correct shape
        else:
            raise ValueError(
                f"Expected 1D or (1, n_features) array, got shape {features.shape}"
            )

        coords = self.umap.transform(features)

        if coords.shape != (1, 2):
            raise RuntimeError(
                f"Unexpected UMAP output shape: {coords.shape}"
            )

        x, y = coords[0]
        return float(x), float(y)

    def transform_batch(self, features: np.ndarray) -> np.ndarray:
        """
        Project multiple molecules into UMAP space.

        Parameters
        ----------
        features : np.ndarray
            Shape (n_samples, n_features)

        Returns
        -------
        np.ndarray
            Shape (n_samples, 2) - UMAP coordinates
        """
        if self.umap is None:
            raise RuntimeError("UMAP model not loaded. Call load() first.")

        if features.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {features.shape}")

        return self.umap.transform(features)


# ------------------------------------------------------------------
# Singleton instance (lazy loaded)
# ------------------------------------------------------------------

_umap_projector: Optional[UMAPProjector] = None


def get_umap_projector() -> UMAPProjector:
    """Get or create the singleton UMAP projector."""
    global _umap_projector
    if _umap_projector is None:
        model_path = Path("backend/ml/models/umap.pkl")
        if model_path.exists():
            _umap_projector = UMAPProjector(str(model_path))
        else:
            raise FileNotFoundError(
                f"UMAP model not found: {model_path}. "
                "Please run training first (see umap_play.ipynb)."
            )
    return _umap_projector


# For backwards compatibility - lazy loading
class _LazyUMAPProjector:
    """Lazy proxy that loads UMAP model on first access."""

    def __getattr__(self, name):
        return getattr(get_umap_projector(), name)


umap_projector = _LazyUMAPProjector()
