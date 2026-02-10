"""
Centralized configuration for the Cardiotoxicity Prediction API.

All configuration values are defined here to avoid hardcoded values
scattered throughout the codebase.
"""

import os
from pathlib import Path
from typing import List

# Base paths
PROJECT_ROOT = Path(__file__).parent
MODELS_DIR = PROJECT_ROOT / "ml" / "models"

# Model file paths
XGB_MODEL_PATH = MODELS_DIR / "xgb_model.pkl"
UMAP_MODEL_PATH = MODELS_DIR / "umap.pkl"
DESC_SCALER_PATH = MODELS_DIR / "desc_scaler.pkl"
UMAP_COORDS_PATH = MODELS_DIR / "umap_coords.csv"

# CORS settings
CORS_ORIGINS: List[str] = os.getenv(
    "CORS_ORIGINS",
    "http://localhost:5173,http://localhost:3000,http://localhost:8000"
).split(",")

# Prediction thresholds
TOXICITY_THRESHOLD = 0.5

# Risk level thresholds for explainability
RISK_THRESHOLDS = {
    "high": 0.8,
    "medium": 0.6,
    "low": 0.4,
}

# Lipinski Rule of 5 thresholds
LIPINSKI_THRESHOLDS = {
    "mw": 500,
    "logp": 5,
    "hbd": 5,
    "hba": 10,
}

# Explainability thresholds
EXPLAINABILITY = {
    "high_logp_threshold": 3.5,
    "high_mw_threshold": 450,
    "many_aromatic_rings_threshold": 3,
    "many_rotatable_bonds_threshold": 8,
}
