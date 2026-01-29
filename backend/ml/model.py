import joblib
import numpy as np
from pathlib import Path
from typing import Optional


class CardiotoxicityModel:
    """
    Wrapper around a pre-trained XGBoost model for cardiotoxicity prediction.
    """

    def __init__(self, model_path: Optional[str] = None):
        """
        Parameters
        ----------
        model_path : str, optional
            Path to serialized model (.pkl)
            If None, model must be loaded later via load()
        """
        self.model = None
        self.scaler = None
        if model_path is not None:
            self.load(model_path)

    def load(self, model_path: str):
        """Load model and scaler from file."""
        path = Path(model_path)
        if not path.exists():
            raise FileNotFoundError(f"Model not found: {model_path}")
        try:
            package = joblib.load(model_path)
            self.model = package["model"]
            self.scaler = package.get("scaler")
        except Exception as e:
            raise RuntimeError(f"Failed to load model from {model_path}") from e

    def predict_proba(self, X: np.ndarray) -> float:
        """
        Predict probability of cardiotoxicity.

        Parameters
        ----------
        X : np.ndarray
            Features with shape (1, n_features) or (n_features,)

        Returns
        -------
        float
            Probability of toxic class (class 1)
        """
        if self.model is None:
            raise RuntimeError("Model not loaded. Call load() first.")

        # Ensure 2D input
        if X.ndim == 1:
            X = X.reshape(1, -1)
        elif X.ndim == 2 and X.shape[0] == 1:
            pass  # Already correct shape
        else:
            raise ValueError(
                f"Expected 1D or (1, n_features) array, got shape {X.shape}"
            )

        proba = self.model.predict_proba(X)
        return float(proba[0, 1])

    def predict_proba_batch(self, X: np.ndarray) -> np.ndarray:
        """
        Predict probabilities for multiple samples.

        Parameters
        ----------
        X : np.ndarray
            Features with shape (n_samples, n_features)

        Returns
        -------
        np.ndarray
            Probabilities of toxic class (class 1), shape (n_samples,)
        """
        if self.model is None:
            raise RuntimeError("Model not loaded. Call load() first.")

        if X.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {X.shape}")

        proba = self.model.predict_proba(X)
        return proba[:, 1]

    def predict(self, X: np.ndarray, threshold: float = 0.5) -> int:
        """
        Predict binary class (0 = non-toxic, 1 = toxic).

        Parameters
        ----------
        X : np.ndarray
            Features with shape (1, n_features) or (n_features,)
        threshold : float
            Classification threshold (default 0.5)

        Returns
        -------
        int
            0 or 1
        """
        prob = self.predict_proba(X)
        return 1 if prob >= threshold else 0


# ------------------------------------------------------------------
# Singleton instance (lazy loaded)
# ------------------------------------------------------------------

_model: Optional[CardiotoxicityModel] = None


def get_model() -> CardiotoxicityModel:
    """Get or create the singleton model."""
    global _model
    if _model is None:
        model_path = Path("backend/ml/models/xgb_model.pkl")
        if model_path.exists():
            _model = CardiotoxicityModel(str(model_path))
        else:
            raise FileNotFoundError(
                f"Model not found: {model_path}. "
                "Please run training first."
            )
    return _model


# For backwards compatibility - lazy loading
class _LazyModel:
    """Lazy proxy that loads model on first access."""

    def __getattr__(self, name):
        return getattr(get_model(), name)


model = _LazyModel()
