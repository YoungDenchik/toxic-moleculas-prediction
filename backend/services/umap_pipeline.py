from typing import Dict

from backend.services.core_utils.standardization import standardize_smiles
from backend.services.core_utils.mol import mol_from_canonical_smiles
from backend.services.core_utils.features import build_feature_vector
from backend.ml.umap_model import umap_projector
from backend.ml.model import model


def project_smiles_to_umap_with_prediction(smiles: str) -> Dict:
    """
    Project a SMILES string to UMAP coordinates and predict toxicity.

    Parameters
    ----------
    smiles : str
        Input SMILES string

    Returns
    -------
    dict
        Contains ok, canonical_smiles, prediction (label + probability), x, y
    """
    canonical, error = standardize_smiles(smiles)

    if error is not None:
        return {
            "ok": False,
            "error": error,
        }

    mol = mol_from_canonical_smiles(canonical)
    if mol is None:
        return {
            "ok": False,
            "error": "Failed to create molecule from SMILES",
        }

    # Build feature vector (returns numpy array of shape (1, n_features))
    features = build_feature_vector(mol)

    # Predict toxicity probability
    prob = model.predict_proba(features)
    label = "toxic" if prob >= 0.5 else "non-toxic"

    # Project to UMAP space
    x, y = umap_projector.transform(features)

    return {
        "ok": True,
        "canonical_smiles": canonical,
        "prediction": {
            "label": label,
            "probability": prob,
        },
        "x": x,
        "y": y,
    }
