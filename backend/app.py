"""
Cardiotoxicity Prediction API

FastAPI application for predicting cardiotoxicity from molecular SMILES.
"""

from typing import List

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from backend.api import predict, chemical_space, properties, analyze, explain
from backend.api.schemas import TrainingUmapPoint
from backend.services.umap_data import load_umap_training_data
from backend.config import CORS_ORIGINS
from backend.exceptions import CardiotoxicityAPIError


app = FastAPI(
    title="Cardiotoxicity Prediction API",
    description="Predict cardiotoxicity of molecules using ML models",
    version="1.0.0",
)


@app.exception_handler(CardiotoxicityAPIError)
async def cardiotoxicity_exception_handler(
    _request: Request, exc: CardiotoxicityAPIError
) -> JSONResponse:
    """Handle custom API exceptions with proper status codes."""
    return JSONResponse(
        status_code=exc.status_code,
        content={"detail": exc.message},
    )

# CORS middleware for frontend integration
app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(predict.router)
app.include_router(chemical_space.router)
app.include_router(properties.router)
app.include_router(analyze.router)
app.include_router(explain.router)


@app.get("/", tags=["Health"])
def root():
    """Health check endpoint."""
    return {"status": "ok", "message": "Cardiotoxicity Prediction API"}


@app.get("/health", tags=["Health"])
def health():
    """Health check endpoint."""
    return {"status": "healthy"}


@app.get("/umap", response_model=List[TrainingUmapPoint], tags=["Chemical Space"])
def get_umap():
    """Get precomputed UMAP coordinates for training data visualization."""
    return load_umap_training_data()


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
