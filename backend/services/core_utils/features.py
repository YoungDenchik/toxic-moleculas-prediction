# # from turtle import pd
# # import numpy as np
# # from backend.services.core_utils.fingerprints import ecfp4_fingerprint
# # from backend.services.core_utils.descriptors import physchem_descriptors




# # def build_feature_vector(mol) -> np.ndarray:
# #     """
# #     Final feature vector used by the ML model.
# #     Must match training pipeline EXACTLY.
# #     """
# #     fp = ecfp4_fingerprint(mol)          # (1024,)
# #     physchem = physchem_descriptors(mol) # (6,)

# #     features = np.concatenate([fp, physchem])
# #     return features


# # # class FeatureBuilder:
# # #     def __init__(
# # #         self,
# # #         maccs_gen,
# # #         scaler,
# # #         descriptor_cols: list[str]
# # #     ):
# # #         self.maccs_gen = maccs_gen
# # #         self.scaler = scaler
# # #         self.descriptor_cols = descriptor_cols

# # #     def transform(self, smiles_list: list[str]) -> pd.DataFrame:
# # #         X_fp = self._maccs(smiles_list)
# # #         X_desc = self._descriptors(smiles_list)
# # #         return pd.concat([X_fp, X_desc], axis=1)


# from turtle import pd
# import numpy as np
# from rdkit.Chem import Descriptors
# from rdkit import Chem

# from skfp.fingerprints import (
#     MACCSFingerprint,
# )

# final_desc_cols = ['MaxAbsEStateIndex',
#  'MinAbsEStateIndex',
#  'MinEStateIndex',
#  'qed',
#  'SPS',
#  'MolWt',
#  'NumRadicalElectrons',
#  'MaxPartialCharge',
#  'MinPartialCharge',
#  'FpDensityMorgan1',
#  'AvgIpc',
#  'BalabanJ',
#  'Ipc',
#  'PEOE_VSA1',
#  'PEOE_VSA10',
#  'PEOE_VSA11',
#  'PEOE_VSA12',
#  'PEOE_VSA13',
#  'PEOE_VSA14',
#  'PEOE_VSA2',
#  'PEOE_VSA3',
#  'PEOE_VSA4',
#  'PEOE_VSA5',
#  'PEOE_VSA6',
#  'PEOE_VSA7',
#  'PEOE_VSA8',
#  'PEOE_VSA9',
#  'SMR_VSA1',
#  'SMR_VSA10',
#  'SMR_VSA2',
#  'SMR_VSA3',
#  'SMR_VSA4',
#  'SMR_VSA5',
#  'SMR_VSA6',
#  'SMR_VSA7',
#  'SMR_VSA9',
#  'SlogP_VSA1',
#  'SlogP_VSA10',
#  'SlogP_VSA11',
#  'SlogP_VSA12',
#  'SlogP_VSA2',
#  'SlogP_VSA3',
#  'SlogP_VSA4',
#  'SlogP_VSA7',
#  'SlogP_VSA8',
#  'TPSA',
#  'EState_VSA1',
#  'EState_VSA11',
#  'EState_VSA2',
#  'EState_VSA3',
#  'EState_VSA4',
#  'EState_VSA5',
#  'EState_VSA6',
#  'EState_VSA7',
#  'EState_VSA8',
#  'EState_VSA9',
#  'VSA_EState2',
#  'VSA_EState3',
#  'VSA_EState4',
#  'VSA_EState5',
#  'VSA_EState7',
#  'VSA_EState8',
#  'VSA_EState9',
#  'FractionCSP3',
#  'NHOHCount',
#  'NumAliphaticCarbocycles',
#  'NumAliphaticHeterocycles',
#  'NumAromaticCarbocycles',
#  'NumAromaticHeterocycles',
#  'NumAromaticRings',
#  'RingCount',
#  'MolLogP',
#  'fr_Al_COO',
#  'fr_Al_OH',
#  'fr_ArN',
#  'fr_Ar_COO',
#  'fr_Ar_NH',
#  'fr_Ar_OH',
#  'fr_C_O',
#  'fr_C_S',
#  'fr_HOCCN',
#  'fr_Imine',
#  'fr_NH1',
#  'fr_NH2',
#  'fr_N_O',
#  'fr_Ndealkylation1',
#  'fr_Ndealkylation2',
#  'fr_SH',
#  'fr_alkyl_carbamate',
#  'fr_allylic_oxid',
#  'fr_amidine',
#  'fr_aniline',
#  'fr_aryl_methyl',
#  'fr_azo',
#  'fr_barbitur',
#  'fr_benzodiazepine',
#  'fr_bicyclic',
#  'fr_dihydropyridine',
#  'fr_epoxide',
#  'fr_ester',
#  'fr_ether',
#  'fr_furan',
#  'fr_guanido',
#  'fr_hdrzine',
#  'fr_hdrzone',
#  'fr_imidazole',
#  'fr_imide',
#  'fr_ketone',
#  'fr_lactam',
#  'fr_lactone',
#  'fr_methoxy',
#  'fr_morpholine',
#  'fr_nitro',
#  'fr_oxazole',
#  'fr_oxime',
#  'fr_para_hydroxylation',
#  'fr_phos_acid',
#  'fr_piperdine',
#  'fr_piperzine',
#  'fr_priamide',
#  'fr_pyridine',
#  'fr_quatN',
#  'fr_sulfide',
#  'fr_sulfonamd',
#  'fr_sulfone',
#  'fr_term_acetylene',
#  'fr_tetrazole',
#  'fr_thiazole',
#  'fr_thiophene',
#  'fr_unbrch_alkane',
#  'fr_urea']

# def create_hybrid(fp_df, desc_df):
#     return pd.concat([fp_df.reset_index(drop=True), 
#                       desc_df.reset_index(drop=True)], axis=1)



# def get_specific_descriptors(smiles_list, descriptor_names):
#     all_funcs = {name: func for name, func in Descriptors._descList}
    
#     selected_funcs = [(name, all_funcs[name]) for name in final_desc_cols if name in all_funcs]
    
#     data = []
#     for smi in smiles_list:
#         mol = Chem.MolFromSmiles(smi)
#         if mol:
#             vals = []
#             for name, f in selected_funcs:
#                 try:
#                     vals.append(f(mol))
#                 except:
#                     vals.append(np.nan)
#         else:
#             vals = [np.nan] * len(selected_funcs)
#         data.append(vals)

#     actual_names = [name for name, _ in selected_funcs]
#     df = pd.DataFrame(data, columns=actual_names)
    
#     return df


# def prepare_test_data(smiles_list, maccs_gen, final_desc_cols, scaler):
#     maccs_features = maccs_gen.transform(smiles_list)
#     X_maccs = pd.DataFrame(maccs_features, columns=[f"MACCS_{i}" for i in range(maccs_features.shape[1])])

#     X_desc_raw = get_specific_descriptors(smiles_list, final_desc_cols)
    
#     X_desc_filtered = X_desc_raw[final_desc_cols]
    
#     X_desc_scaled_np = scaler.transform(X_desc_filtered)
#     X_desc_scaled = pd.DataFrame(X_desc_scaled_np, columns=final_desc_cols)
    
#     X_hybrid = create_hybrid(X_maccs, X_desc_scaled)
    
#     return X_hybrid

# maccs_gen = MACCSFingerprint(n_jobs=-1)


import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

from skfp.fingerprints import MACCSFingerprint


# maccs_gen = MACCSFingerprint(n_jobs=-1)

# scaler_path = joblib.load("backend/ml/models/xgb_model.pkl")
# scaler = scaler_path["scaler"]
# model = joblib.load("tox_model.joblib")
# umap_projector = joblib.load("umap.joblib")



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


def get_maccs_for_mol(
    mol: Chem.Mol,
    maccs_gen: MACCSFingerprint,
) -> pd.DataFrame:
    """
    Returns MACCS fingerprint for ONE molecule as DataFrame (1 × 166)
    """
    fp = maccs_gen.transform([mol])  # IMPORTANT: list!
    cols = [f"MACCS_{i}" for i in range(fp.shape[1])]
    return pd.DataFrame(fp, columns=cols)


def create_hybrid_features(
    X_fp: pd.DataFrame,
    X_desc: pd.DataFrame,
) -> pd.DataFrame:
    return pd.concat([X_fp.reset_index(drop=True),
                      X_desc.reset_index(drop=True)], axis=1)


def build_feature_vector(
    mol: Chem.Mol,
    # *,
    # maccs_gen: MACCSFingerprint,
    # scaler,
    # final_desc_cols: list[str],
) -> pd.DataFrame:
    """
    Build full feature vector for ONE molecule.
    Output shape: (1, N_features)
    """

    final_desc_cols = ['MaxAbsEStateIndex',
    'MinAbsEStateIndex',
    'MinEStateIndex',
    'qed',
    'SPS',
    'MolWt',
    'NumRadicalElectrons',
    'MaxPartialCharge',
    'MinPartialCharge',
    'FpDensityMorgan1',
    'AvgIpc',
    'BalabanJ',
    'Ipc',
    'PEOE_VSA1',
    'PEOE_VSA10',
    'PEOE_VSA11',
    'PEOE_VSA12',
    'PEOE_VSA13',
    'PEOE_VSA14',
    'PEOE_VSA2',
    'PEOE_VSA3',
    'PEOE_VSA4',
    'PEOE_VSA5',
    'PEOE_VSA6',
    'PEOE_VSA7',
    'PEOE_VSA8',
    'PEOE_VSA9',
    'SMR_VSA1',
    'SMR_VSA10',
    'SMR_VSA2',
    'SMR_VSA3',
    'SMR_VSA4',
    'SMR_VSA5',
    'SMR_VSA6',
    'SMR_VSA7',
    'SMR_VSA9',
    'SlogP_VSA1',
    'SlogP_VSA10',
    'SlogP_VSA11',
    'SlogP_VSA12',
    'SlogP_VSA2',
    'SlogP_VSA3',
    'SlogP_VSA4',
    'SlogP_VSA7',
    'SlogP_VSA8',
    'TPSA',
    'EState_VSA1',
    'EState_VSA11',
    'EState_VSA2',
    'EState_VSA3',
    'EState_VSA4',
    'EState_VSA5',
    'EState_VSA6',
    'EState_VSA7',
    'EState_VSA8',
    'EState_VSA9',
    'VSA_EState2',
    'VSA_EState3',
    'VSA_EState4',
    'VSA_EState5',
    'VSA_EState7',
    'VSA_EState8',
    'VSA_EState9',
    'FractionCSP3',
    'NHOHCount',
    'NumAliphaticCarbocycles',
    'NumAliphaticHeterocycles',
    'NumAromaticCarbocycles',
    'NumAromaticHeterocycles',
    'NumAromaticRings',
    'RingCount',
    'MolLogP',
    'fr_Al_COO',
    'fr_Al_OH',
    'fr_ArN',
    'fr_Ar_COO',
    'fr_Ar_NH',
    'fr_Ar_OH',
    'fr_C_O',
    'fr_C_S',
    'fr_HOCCN',
    'fr_Imine',
    'fr_NH1',
    'fr_NH2',
    'fr_N_O',
    'fr_Ndealkylation1',
    'fr_Ndealkylation2',
    'fr_SH',
    'fr_alkyl_carbamate',
    'fr_allylic_oxid',
    'fr_amidine',
    'fr_aniline',
    'fr_aryl_methyl',
    'fr_azo',
    'fr_barbitur',
    'fr_benzodiazepine',
    'fr_bicyclic',
    'fr_dihydropyridine',
    'fr_epoxide',
    'fr_ester',
    'fr_ether',
    'fr_furan',
    'fr_guanido',
    'fr_hdrzine',
    'fr_hdrzone',
    'fr_imidazole',
    'fr_imide',
    'fr_ketone',
    'fr_lactam',
    'fr_lactone',
    'fr_methoxy',
    'fr_morpholine',
    'fr_nitro',
    'fr_oxazole',
    'fr_oxime',
    'fr_para_hydroxylation',
    'fr_phos_acid',
    'fr_piperdine',
    'fr_piperzine',
    'fr_priamide',
    'fr_pyridine',
    'fr_quatN',
    'fr_sulfide',
    'fr_sulfonamd',
    'fr_sulfone',
    'fr_term_acetylene',
    'fr_tetrazole',
    'fr_thiazole',
    'fr_thiophene',
    'fr_unbrch_alkane',
    'fr_urea']

    maccs_gen = MACCSFingerprint(n_jobs=-1)

    scaler_path = joblib.load("backend/ml/models/xgb_model.pkl")
    scaler = scaler_path["scaler"]

    # 1. MACCS
    X_maccs = get_maccs_for_mol(mol, maccs_gen)

    # 2. Descriptors (raw)
    X_desc_raw = get_specific_descriptors_for_mol(mol, final_desc_cols)

    # 3. Scale descriptors (IMPORTANT: same scaler as train)
    X_desc_scaled_np = scaler.transform(X_desc_raw)
    X_desc_scaled = pd.DataFrame(
        X_desc_scaled_np,
        columns=final_desc_cols,
    )

    # 4. Hybrid
    X = create_hybrid_features(X_maccs, X_desc_scaled)

    return X
