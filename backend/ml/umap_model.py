import joblib
import numpy as np
from typing import Tuple


class UMAPProjector:
    """
    Wrapper around a pre-trained UMAP model.
    Responsible ONLY for projection (transform), never fit.
    """

    def __init__(self, model_path: str):
        """
        Parameters
        ----------
        model_path : str
            Path to serialized UMAP model (.pkl)
        """
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
            Shape (n_features,)

        Returns
        -------
        (x, y) : Tuple[float, float]
            2D UMAP coordinates
        """
        if features.ndim != 1:
            raise ValueError(
                f"Expected 1D feature vector, got shape {features.shape}"
            )

        # UMAP expects 2D input: (n_samples, n_features)
        features_2d = features.reshape(1, -1)

        coords = self.umap.transform(features_2d)

        if coords.shape != (1, 2):
            raise RuntimeError(
                f"Unexpected UMAP output shape: {coords.shape}"
            )

        x, y = coords[0]
        return float(x), float(y)


# ------------------------------------------------------------------
# Singleton instance (loaded once at startup)
# ------------------------------------------------------------------

umap_projector = UMAPProjector(
    model_path="backend/ml/models/umap.pkl"
)
