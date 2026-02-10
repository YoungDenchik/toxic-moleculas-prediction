"""
Molecule Service

Centralized service for molecule processing, feature extraction, and prediction.
Eliminates code duplication between prediction_pipeline.py and umap_pipeline.py.
"""

from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np

from backend.services.core_utils.standardization import standardize_smiles
from backend.services.core_utils.mol import mol_from_canonical_smiles
from backend.services.core_utils.features import build_feature_vector
from backend.ml.model import model
from backend.ml.umap_model import umap_projector
from backend.config import TOXICITY_THRESHOLD


@dataclass
class MoleculeData:
    """Container for processed molecule data."""
    canonical_smiles: str
    mol: object  # RDKit Mol object
    features: np.ndarray


@dataclass
class PredictionResult:
    """Container for prediction results."""
    label: str
    probability: float


def process_smiles(smiles: str) -> Tuple[Optional[MoleculeData], Optional[str]]:
    """
    Process a SMILES string: standardize, create mol object, and extract features.

    Parameters
    ----------
    smiles : str
        Input SMILES string

    Returns
    -------
    Tuple[Optional[MoleculeData], Optional[str]]
        (MoleculeData, None) on success, (None, error_message) on failure
    """
    canonical, error = standardize_smiles(smiles)
    if error is not None:
        return None, error

    mol = mol_from_canonical_smiles(canonical)
    if mol is None:
        return None, "Failed to create molecule from SMILES"

    features = build_feature_vector(mol)

    return MoleculeData(
        canonical_smiles=canonical,
        mol=mol,
        features=features
    ), None


def predict_toxicity(features: np.ndarray) -> PredictionResult:
    """
    Predict toxicity from feature vector.

    Parameters
    ----------
    features : np.ndarray
        Feature vector of shape (1, n_features)

    Returns
    -------
    PredictionResult
        Prediction label and probability
    """
    prob = model.predict_proba(features)
    label = "toxic" if prob >= TOXICITY_THRESHOLD else "non-toxic"
    return PredictionResult(label=label, probability=prob)


def project_to_umap(features: np.ndarray) -> Tuple[float, float]:
    """
    Project features to UMAP coordinates.

    Parameters
    ----------
    features : np.ndarray
        Feature vector of shape (1, n_features)

    Returns
    -------
    Tuple[float, float]
        (x, y) UMAP coordinates
    """
    return umap_projector.transform(features)
