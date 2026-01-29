from typing import Dict

from backend.services.core_utils.standardization import standardize_smiles
from backend.services.core_utils.mol import mol_from_canonical_smiles
from backend.services.core_utils.physchem import compute_physchem, check_lipinski


def analyze_lipinski(smiles: str) -> Dict:
    canonical, error = standardize_smiles(smiles)

    if error is not None:
        return {
            "ok": False,
            "error": error,
        }

    mol = mol_from_canonical_smiles(canonical)

    props = compute_physchem(mol)
    lipinski = check_lipinski(props)

    return {
        "ok": True,
        "canonical_smiles": canonical,
        "properties": props,
        "checks": lipinski["checks"],
        "lipinski_pass": lipinski["pass"],
    }

