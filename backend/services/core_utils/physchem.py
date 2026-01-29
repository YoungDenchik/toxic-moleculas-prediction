from rdkit.Chem import Descriptors


def compute_physchem(mol) -> dict:
    """
    Compute physicochemical properties required for Lipinski Rule of 5.
    """
    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
    }


def check_lipinski(props: dict) -> dict:
    """
    Check Lipinski Rule of 5.
    """
    checks = {
        "molecular_weight": props["molecular_weight"] <= 500,
        "logp": props["logp"] <= 5,
        "hbd": props["hbd"] <= 5,
        "hba": props["hba"] <= 10,
    }

    return {
        "checks": checks,
        "pass": all(checks.values()),
    }
