from typing import Dict

from backend.services.core_utils.standardization import standardize_smiles
from backend.services.core_utils.mol import mol_from_canonical_smiles
from backend.services.core_utils.features import build_feature_vector
from backend.ml.model import model


def predict_from_smiles(smiles: str) -> Dict:
    canonical, error = standardize_smiles(smiles)

    if error is not None:
        return {
            "ok": False,
            "error": error,
        }

    mol = mol_from_canonical_smiles(canonical)

    features = build_feature_vector(mol)

    prob = model.predict_proba(features)
    label = "toxic" if prob >= 0.5 else "non-toxic"

    return {
        "ok": True,
        "canonical_smiles": canonical,
        "prediction": {
            "label": label,
            "probability": prob,
        },
    }
