from fastapi import APIRouter, HTTPException

from backend.api.schemas import (
    PhysChemPropertiesResponse,
    PropertyResult,
)
from backend.api.schemas import PredictionRequest
from backend.services.lipinski_pipeline import analyze_lipinski

router = APIRouter(
    prefix="/properties",
    tags=["Physicochemical Properties"],
)


@router.post(
    "/lipinski",
    response_model=PhysChemPropertiesResponse,
    summary="Calculate physicochemical properties and Lipinski Rule of 5",
)
def lipinski_properties(request: PredictionRequest):
    """
    SMILES -> physicochemical properties -> Lipinski Rule of 5
    """
    result = analyze_lipinski(request.smiles)

    if not result["ok"]:
        raise HTTPException(
            status_code=400,
            detail=result["error"],
        )

    props = result["properties"]
    checks = result["checks"]

    return PhysChemPropertiesResponse(
        canonical_smiles=result["canonical_smiles"],
        molecular_weight=PropertyResult(
            value=props["molecular_weight"],
            pass_rule=checks["molecular_weight"],
        ),
        logp=PropertyResult(
            value=props["logp"],
            pass_rule=checks["logp"],
        ),
        hbd=PropertyResult(
            value=props["hbd"],
            pass_rule=checks["hbd"],
        ),
        hba=PropertyResult(
            value=props["hba"],
            pass_rule=checks["hba"],
        ),
        lipinski_pass=result["lipinski_pass"],
    )
