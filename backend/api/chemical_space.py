from fastapi import APIRouter, HTTPException

from backend.api.schemas import (
    UMAPWithPredictionResponse,
    UMAPPoint,
)
from backend.api.schemas import PredictionOutput
from backend.api.schemas import PredictionRequest
from backend.services.umap_pipeline import project_smiles_to_umap_with_prediction

router = APIRouter(
    prefix="/chemical-space",
    tags=["Chemical Space"],
)


@router.post(
    "/umap",
    response_model=UMAPWithPredictionResponse,
    summary="Project molecule into UMAP space with prediction",
)
def project_umap(request: PredictionRequest):
    """
    SMILES -> features -> prediction -> UMAP.transform
    """
    result = project_smiles_to_umap_with_prediction(request.smiles)

    if not result["ok"]:
        raise HTTPException(
            status_code=400,
            detail=result["error"],
        )

    return UMAPWithPredictionResponse(
        canonical_smiles=result["canonical_smiles"],
        prediction=result["prediction"],
        point=UMAPPoint(
            x=result["x"],
            y=result["y"],
        ),
    )
