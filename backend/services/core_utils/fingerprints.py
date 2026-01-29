import numpy as np
import pandas
from rdkit.Chem import AllChem, MACCSkeys
import pandas as pd

def ecfp4_fingerprint(mol, n_bits: int = 1024) -> np.ndarray:
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol,
        radius=2,
        nBits=n_bits,
    )
    return np.asarray(fp, dtype=np.float32)


def maccs_fingerprint(mol) -> np.ndarray:
    fp = MACCSkeys.GenMACCSKeys(mol)
    return np.asarray(fp, dtype=np.float32)


class MACCSFeature:
    def __init__(self, generator):
        self.generator = generator

    def transform(self, smiles_list):
        X = self.generator.transform(smiles_list)
        return pd.DataFrame(
            X,
            columns=[f"MACCS_{i}" for i in range(X.shape[1])]
        )
