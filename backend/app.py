"""
Cardiotoxicity Prediction API

FastAPI application for predicting cardiotoxicity from molecular SMILES.
"""

from typing import List

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.api import predict, chemical_space, properties, analyze, explain
from backend.api.schemas import TrainingUmapPoint
from backend.services.umap_data import load_umap_training_data


app = FastAPI(
    title="Cardiotoxicity Prediction API",
    description="Predict cardiotoxicity of molecules using ML models",
    version="1.0.0",
)

# CORS middleware for frontend integration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
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
