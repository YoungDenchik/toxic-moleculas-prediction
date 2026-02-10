"""
Prediction Pipeline

Predicts cardiotoxicity from SMILES input.
"""

from typing import Dict

from backend.services.molecule_service import process_smiles, predict_toxicity


def predict_from_smiles(smiles: str) -> Dict:
    """
    Predict cardiotoxicity from a SMILES string.

    Parameters
    ----------
    smiles : str
        Input SMILES string

    Returns
    -------
    Dict
        Contains ok, canonical_smiles, prediction (label + probability)
    """
    mol_data, error = process_smiles(smiles)

    if error is not None:
        return {
            "ok": False,
            "error": error,
        }

    prediction = predict_toxicity(mol_data.features)

    return {
        "ok": True,
        "canonical_smiles": mol_data.canonical_smiles,
        "prediction": {
            "label": prediction.label,
            "probability": prediction.probability,
        },
    }
