from fastapi import APIRouter, HTTPException

from backend.api.schemas import (
    PredictionRequest,
    PredictionResponse,
)
from backend.services.prediction_pipeline import predict_from_smiles

router = APIRouter(
    prefix="/predict",
    tags=["Prediction"],
)


@router.post(
    "",
    response_model=PredictionResponse,
    summary="Predict cardiotoxicity from SMILES",
)
def predict(request: PredictionRequest):
    """
    Prediction pipeline:
    raw SMILES -> standardization -> features -> ML model -> prediction
    """
    result = predict_from_smiles(request.smiles)

    if not result["ok"]:
        raise HTTPException(
            status_code=400,
            detail=result["error"],
        )

    return PredictionResponse(
        canonical_smiles=result["canonical_smiles"],
        prediction=result["prediction"],
    )
