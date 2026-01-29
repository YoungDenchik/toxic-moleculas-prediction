from typing import Optional, Tuple
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def standardize_smiles(
    smiles: str,
    remove_stereo: bool = False,
) -> Tuple[Optional[str], Optional[str]]:

    if not isinstance(smiles, str) or smiles.strip() == "":
        return None, "empty_smiles"

    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception:
        return None, "rdkit_parse_error"

    if mol is None:
        return None, "invalid_smiles"

    try:
        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)
    except Exception:
        return None, "fragment_chooser_error"

    try:
        taut_enum = rdMolStandardize.TautomerEnumerator()
        mol = taut_enum.Canonicalize(mol)
    except Exception:
        return None, "tautomer_error"

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None, "sanitize_error"

    try:
        canonical = Chem.MolToSmiles(
            mol,
            canonical=True,
            isomericSmiles=not remove_stereo
        )
    except Exception:
        return None, "smiles_generation_error"

    return canonical, None
