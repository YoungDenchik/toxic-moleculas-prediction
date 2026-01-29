import joblib
import numpy as np
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors

from skfp.fingerprints import MACCSFingerprint

from backend.services.core_utils.descriptors import final_desc_cols


# Load artifacts once at module level
_ARTIFACTS_DIR = Path("backend/ml/models")
_maccs_gen = MACCSFingerprint(n_jobs=-1)
_scaler = None


def _get_scaler():
    """Lazily load scaler to avoid errors if model file doesn't exist yet."""
    global _scaler
    if _scaler is None:
        scaler_path = _ARTIFACTS_DIR / "xgb_model.pkl"
        if scaler_path.exists():
            package = joblib.load(scaler_path)
            _scaler = package["scaler"]
        else:
            raise FileNotFoundError(
                f"Model file not found: {scaler_path}. "
                "Please run training first."
            )
    return _scaler


def get_specific_descriptors_for_mol(
    mol: Chem.Mol,
    descriptor_names: list[str],
) -> pd.DataFrame:
    """
    Compute selected RDKit descriptors for ONE molecule.
    Returns DataFrame with shape (1, n_descriptors)
    """
    all_funcs = dict(Descriptors._descList)

    values = {}
    for name in descriptor_names:
        func = all_funcs.get(name)
        if func is None:
            values[name] = np.nan
            continue
        try:
            values[name] = func(mol)
        except Exception:
            values[name] = np.nan

    return pd.DataFrame([values])


def get_maccs_for_mol(mol: Chem.Mol) -> pd.DataFrame:
    """
    Returns MACCS fingerprint for ONE molecule as DataFrame (1 × 167)
    """
    fp = _maccs_gen.transform([mol])
    cols = [f"MACCS_{i}" for i in range(fp.shape[1])]
    return pd.DataFrame(fp, columns=cols)


def create_hybrid_features(
    X_fp: pd.DataFrame,
    X_desc: pd.DataFrame,
) -> pd.DataFrame:
    return pd.concat([X_fp.reset_index(drop=True),
                      X_desc.reset_index(drop=True)], axis=1)


def build_feature_vector(mol: Chem.Mol) -> np.ndarray:
    """
    Build full feature vector for ONE molecule.
    Output shape: (1, N_features) as numpy array

    Features: MACCS fingerprints (167) + scaled descriptors (130)
    """
    scaler = _get_scaler()

    # 1. MACCS fingerprints
    X_maccs = get_maccs_for_mol(mol)

    # 2. Descriptors (raw)
    X_desc_raw = get_specific_descriptors_for_mol(mol, final_desc_cols)

    # 3. Scale descriptors (same scaler as training)
    X_desc_scaled_np = scaler.transform(X_desc_raw)
    X_desc_scaled = pd.DataFrame(
        X_desc_scaled_np,
        columns=final_desc_cols,
    )

    # 4. Create hybrid features
    X = create_hybrid_features(X_maccs, X_desc_scaled)

    return X.values


def build_features_for_smiles_list(
    smiles_list: list[str],
    scaler=None,
    fit_scaler: bool = False
) -> pd.DataFrame:
    """
    Build feature matrix for a list of SMILES strings.
    Used for training.

    Parameters
    ----------
    smiles_list : list[str]
        List of SMILES strings
    scaler : StandardScaler, optional
        Scaler to use. If None and fit_scaler=False, uses pre-trained scaler.
    fit_scaler : bool
        If True, fit the scaler on this data (for training)

    Returns
    -------
    pd.DataFrame
        Feature matrix with shape (n_samples, n_features)
    """
    from sklearn.preprocessing import StandardScaler

    if scaler is None and not fit_scaler:
        scaler = _get_scaler()
    elif scaler is None and fit_scaler:
        scaler = StandardScaler()

    # MACCS fingerprints
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    X_maccs = _maccs_gen.transform(mols)
    X_maccs = pd.DataFrame(
        X_maccs,
        columns=[f"MACCS_{i}" for i in range(X_maccs.shape[1])]
    )

    # Descriptors
    all_funcs = dict(Descriptors._descList)
    data = []
    for mol in mols:
        if mol is None:
            data.append([np.nan] * len(final_desc_cols))
            continue
        row = []
        for name in final_desc_cols:
            func = all_funcs.get(name)
            if func is None:
                row.append(np.nan)
                continue
            try:
                row.append(func(mol))
            except Exception:
                row.append(np.nan)
        data.append(row)

    X_desc_raw = pd.DataFrame(data, columns=final_desc_cols)

    # Scale
    if fit_scaler:
        X_desc_scaled_np = scaler.fit_transform(X_desc_raw)
    else:
        X_desc_scaled_np = scaler.transform(X_desc_raw)

    X_desc_scaled = pd.DataFrame(X_desc_scaled_np, columns=final_desc_cols)

    # Combine
    X = pd.concat([X_maccs, X_desc_scaled], axis=1)

    return X, scaler
