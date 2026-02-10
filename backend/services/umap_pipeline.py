"""
UMAP Pipeline

Projects molecules to UMAP space and predicts toxicity.
"""

from typing import Dict

from backend.services.molecule_service import (
    process_smiles,
    predict_toxicity,
    project_to_umap,
)


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
    mol_data, error = process_smiles(smiles)

    if error is not None:
        return {
            "ok": False,
            "error": error,
        }

    prediction = predict_toxicity(mol_data.features)
    x, y = project_to_umap(mol_data.features)

    return {
        "ok": True,
        "canonical_smiles": mol_data.canonical_smiles,
        "prediction": {
            "label": prediction.label,
            "probability": prediction.probability,
        },
        "x": x,
        "y": y,
    }
