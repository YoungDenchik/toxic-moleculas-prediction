from rdkit import Chem


def mol_from_canonical_smiles(canonical_smiles: str):
    mol = Chem.MolFromSmiles(canonical_smiles)
    if mol is None:
        raise ValueError("Failed to create RDKit Mol from canonical SMILES")
    return mol
