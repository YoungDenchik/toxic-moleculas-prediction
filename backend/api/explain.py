from fastapi import APIRouter, HTTPException
from backend.api.schemas import PredictionRequest
from backend.api.schemas import ExplainabilityResponse
from backend.services.explainability import explain_toxicity
from backend.services import summarize_with_llm
from backend.services.core_utils.standardization import standardize_smiles
from backend.services.core_utils.mol import mol_from_canonical_smiles
from backend.services.prediction_pipeline import predict_from_smiles

router = APIRouter(
    prefix="/explain",
    tags=["Explainability"],
)


@router.post(
    "/toxicity",
    response_model=ExplainabilityResponse,
    summary="Explain why the molecule might be toxic",
)
def explain(request: PredictionRequest):

    pred = predict_from_smiles(request.smiles)
    if not pred["ok"]:
        raise HTTPException(status_code=400, detail=pred["error"])

    if pred["prediction"]["label"] != "toxic":
        raise HTTPException(
            status_code=400,
            detail="Explainability is only available for toxic predictions",
        )

    canonical = pred["canonical_smiles"]
    mol = mol_from_canonical_smiles(canonical)

    explanation = explain_toxicity(
        mol=mol,
        prediction_probability=pred["prediction"]["probability"],
    )

    # llm_summary = summarize_with_llm(explanation["factors"])

    return ExplainabilityResponse(
        probability=explanation["probability"],
        factors=explanation["factors"],
        # llm_summary=llm_summary,
    )
