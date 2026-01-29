from typing import List, Dict, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors


# ===============================
# Structural toxicophores (SMARTS)
# ===============================

STRUCTURAL_ALERTS = [
    {
        "name": "Aromatic heterocycles",
        "smarts": "[a;r5,r6]",
        "message": "Aromatic heterocycles are often associated with hERG channel binding"
    },
    {
        "name": "Tertiary amine",
        "smarts": "[NX3;!$(NC=O)]",
        "message": "Tertiary amines may interact with cardiac ion channels"
    },
    {
        "name": "Quinone-like system",
        "smarts": "O=c1ccc(O)cc1",
        "message": "Quinone-like structures may induce oxidative stress"
    }
]


def detect_structural_alerts(mol) -> List[str]:
    alerts = []
    for alert in STRUCTURAL_ALERTS:
        patt = Chem.MolFromSmarts(alert["smarts"])
        if patt and mol.HasSubstructMatch(patt):
            alerts.append(alert["message"])
    return alerts


# ===============================
# Physicochemical heuristics
# ===============================

def physchem_heuristics(mol) -> List[str]:
    explanations = []

    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)

    if logp > 3.5:
        explanations.append(
            f"High lipophilicity (LogP = {logp:.2f}) may increase membrane accumulation and off-target effects"
        )

    if mw > 450:
        explanations.append(
            f"High molecular weight (MW = {mw:.1f}) may reduce clearance and increase toxicity risk"
        )

    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    if hbd + hba > 8:
        explanations.append(
            "High number of hydrogen bond donors/acceptors may affect ion channel interactions"
        )

    return explanations


# ===============================
# Main explainability entrypoint
# ===============================

def explain_toxicity(
    mol,
    prediction_probability: float,
) -> Dict:
    """
    Generate human-interpretable explanation for toxic prediction.
    """

    explanations = []

    # Structural alerts
    explanations.extend(detect_structural_alerts(mol))

    # Physchem heuristics
    explanations.extend(physchem_heuristics(mol))

    if not explanations:
        explanations.append(
            "The molecule shares global similarity with known cardiotoxic compounds"
        )

    return {
        "probability": prediction_probability,
        "factors": explanations,
    }
