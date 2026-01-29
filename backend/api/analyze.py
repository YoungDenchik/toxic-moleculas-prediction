from fastapi import APIRouter, HTTPException

from backend.api.schemas import ExplainabilityResponse, PredictionRequest
from backend.api.schemas import AnalyzeResponse
from backend.api.schemas import PropertyResult
from backend.api.schemas import UMAPPoint

from backend.services.core_utils.mol import mol_from_canonical_smiles
from backend.services.explainability import explain_toxicity
from backend.services.prediction_pipeline import predict_from_smiles
from backend.services.umap_pipeline import project_smiles_to_umap_with_prediction
from backend.services.lipinski_pipeline import analyze_lipinski

router = APIRouter(
    prefix="/analyze",
    tags=["Analyze"],
)

@router.post(
    "",
    response_model=AnalyzeResponse,
    summary="Full molecule analysis (prediction + UMAP + Lipinski + explainability)",
)
def analyze(request: PredictionRequest):
    """
    Full analysis pipeline for a molecule.
    """

    # # 1️ Prediction
    # pred = predict_from_smiles(request.smiles)
    # if not pred["ok"]:
    #     raise HTTPException(status_code=400, detail=pred["error"])

    # canonical_smiles = pred["canonical_smiles"]
    # prediction = pred["prediction"]

    # 2️ UMAP
    umap = project_smiles_to_umap_with_prediction(request.smiles)
    if not umap["ok"]:
        raise HTTPException(status_code=400, detail=umap["error"])
    
    canonical_smiles = umap["canonical_smiles"]
    prediction = umap["prediction"]

    # 3️ Lipinski / physchem
    lip = analyze_lipinski(canonical_smiles)
    if not lip["ok"]:
        raise HTTPException(status_code=400, detail=lip["error"])

    props = lip["properties"]
    checks = lip["checks"]

    # 4️ Explainability (ONLY if toxic)
    explainability = None
    if prediction["label"] == "toxic":
        mol = mol_from_canonical_smiles(canonical_smiles)

        explanation = explain_toxicity(
            mol=mol,
            prediction_probability=prediction["probability"],
        )

        # llm_summary = summarize_with_llm(explanation["factors"])

        explainability = ExplainabilityResponse(
            probability=prediction["probability"],
            factors=explanation["factors"],
            # llm_summary=llm_summary,
        )

    return AnalyzeResponse(
        canonical_smiles=canonical_smiles,

        prediction=prediction,

        chemical_space=UMAPPoint(
            x=umap["x"],
            y=umap["y"],
        ),

        properties={
            "molecular_weight": PropertyResult(
                value=props["molecular_weight"],
                pass_rule=checks["molecular_weight"],
            ),
            "logp": PropertyResult(
                value=props["logp"],
                pass_rule=checks["logp"],
            ),
            "hbd": PropertyResult(
                value=props["hbd"],
                pass_rule=checks["hbd"],
            ),
            "hba": PropertyResult(
                value=props["hba"],
                pass_rule=checks["hba"],
            ),
        },

        lipinski_pass=lip["lipinski_pass"],

        explainability=explainability,
    )
