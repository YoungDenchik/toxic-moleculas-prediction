from typing import Dict

from backend.services.core_utils.standardization import standardize_smiles
from backend.services.core_utils.mol import mol_from_canonical_smiles
from backend.services.core_utils.features import build_feature_vector
from backend.ml.umap_model import umap_projector
from backend.ml.model import model

from skfp.fingerprints import MACCSFingerprint


maccs_gen = MACCSFingerprint(n_jobs=-1)

# scaler = joblib.load("descriptor_scaler.joblib")
# model = joblib.load("tox_model.joblib")
# umap_projector = joblib.load("umap.joblib")



def project_smiles_to_umap_with_prediction(smiles: str) -> Dict:
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


# def project_smiles_to_umap_with_prediction(smiles: str) -> Dict:
#     canonical, error = standardize_smiles(smiles)
#     if error is not None:
#         return {"ok": False, "error": error}

#     mol = mol_from_canonical_smiles(canonical)

#     features = build_feature_vector(
#         mol,
#         maccs_gen=maccs_gen,
#         scaler=scaler,
#         final_desc_cols=final_desc_cols,
#     )

#     prob = float(model.predict_proba(features)[0, 1])
#     label = "toxic" if prob >= 0.5 else "non-toxic"

#     x, y = umap_projector.transform(features)[0]

#     return {
#         "ok": True,
#         "canonical_smiles": canonical,
#         "prediction": {
#             "label": label,
#             "probability": prob,
#         },
#         "x": float(x),
#         "y": float(y),
#     }
