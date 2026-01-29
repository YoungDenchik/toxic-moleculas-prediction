"""
Cardiotoxicity Prediction API

FastAPI application for predicting cardiotoxicity from molecular SMILES.
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.api import predict, chemical_space, properties, analyze, explain


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


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
