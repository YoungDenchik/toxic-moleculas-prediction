from pydantic import BaseModel, Field


# ===== Request =====

class PredictionRequest(BaseModel):
    smiles: str = Field(
        ...,
        description="Input molecule as SMILES",
        example="CCOc1ccc2nc(S(=O)(=O)N)sc2c1",
    )


# ===== Response =====

class PredictionOutput(BaseModel):
    label: str = Field(
        ...,
        description="Predicted class label",
        example="toxic",
    )
    probability: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Predicted probability",
        example=0.76,
    )


class PredictionResponse(BaseModel):
    canonical_smiles: str = Field(
        ...,
        description="Canonical SMILES after standardization",
        example="CCOc1ccc2nc(S(=O)(=O)N)sc2c1",
    )
    prediction: PredictionOutput


class UMAPPoint(BaseModel):
    x: float
    y: float


class UMAPWithPredictionResponse(BaseModel):
    method: str = Field(default="UMAP")
    canonical_smiles: str
    prediction: PredictionOutput
    point: UMAPPoint


from pydantic import BaseModel, Field


class PropertyResult(BaseModel):
    value: float
    pass_rule: bool


class PhysChemPropertiesResponse(BaseModel):
    canonical_smiles: str

    molecular_weight: PropertyResult = Field(..., description="MW ≤ 500")
    logp: PropertyResult = Field(..., description="LogP ≤ 5")
    hbd: PropertyResult = Field(..., description="H-bond donors ≤ 5")
    hba: PropertyResult = Field(..., description="H-bond acceptors ≤ 10")

    lipinski_pass: bool

from pydantic import BaseModel, Field
from typing import List, Optional


class ExplainabilityResponse(BaseModel):
    probability: float = Field(..., ge=0.0, le=1.0)
    factors: List[str]
    # llm_summary: Optional[str] = None


class AnalyzeResponse(BaseModel):
    canonical_smiles: str

    prediction: PredictionOutput

    chemical_space: UMAPPoint

    properties: dict[str, PropertyResult]

    lipinski_pass: bool

    explainability: Optional[ExplainabilityResponse] = None



