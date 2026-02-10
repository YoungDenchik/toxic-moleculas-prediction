# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Full-stack ML web application for predicting hERG cardiotoxicity of drug molecules. Users draw or input molecular structures (SMILES), and the app predicts toxicity using an XGBoost model, visualizes the molecule in chemical space via UMAP, and provides structural explanations.

## Build and Run Commands

### Development (Docker with hot reload)
```bash
docker-compose -f docker-compose.dev.yml up
# Backend: http://localhost:8000
# Frontend: http://localhost:5173
```

### Production (Docker)
```bash
docker-compose up
# Frontend: http://localhost:3000 (nginx)
# Backend: http://localhost:8000
```

### Local Development (no Docker)
```bash
# Backend
cd backend
pip install -r requirements.txt
uvicorn app:app --reload --host 0.0.0.0 --port 8000

# Frontend
cd frontend
npm install
npm run dev
```

### Testing
```bash
cd backend
pytest                    # Run all tests
pytest tests/test_predict.py  # Run specific test file
pytest -v                 # Verbose output
```

### Frontend Commands
```bash
npm run dev      # Start Vite dev server
npm run build    # TypeScript compile + production build
npm run lint     # ESLint check
npm run preview  # Preview production build
```

## Architecture

### Backend (`backend/`)
- **Framework:** FastAPI with Uvicorn
- **Entry point:** `app.py`
- **Configuration:** `config.py` - Centralized configuration (paths, thresholds, CORS)
- **Exceptions:** `exceptions.py` - Custom API exceptions

**API Routers (`api/`):**
- `predict.py` - `/predict` - SMILES → toxic/non-toxic prediction
- `chemical_space.py` - `/chemical-space/umap` - SMILES → 2D UMAP coordinates
- `properties.py` - `/properties/lipinski` - SMILES → Lipinski Rule of 5 properties
- `analyze.py` - `/analyze` - Full analysis combining all endpoints
- `explain.py` - `/explain/toxicity` - Structural alerts and risk factors

**Service Layer (`services/`):**
- `molecule_service.py` - Centralized molecule processing (standardization, features, prediction)
- `prediction_pipeline.py` - End-to-end prediction workflow
- `umap_pipeline.py` - UMAP projection with predictions
- `lipinski_pipeline.py` - Physicochemical property analysis
- `explainability.py` - Structural alert detection
- `umap_data.py` - UMAP training data loading

**Core Utilities (`services/core_utils/`):**
- `standardization.py` - SMILES canonicalization
- `mol.py` - RDKit molecule creation
- `features.py` - Feature vector building (167 MACCS + 130 descriptors = 297 total)
- `descriptors.py` - Descriptor column configuration

**ML Models (`ml/`):**
- `model.py` - XGBoost classifier wrapper (lazy loading singleton)
- `umap_model.py` - UMAP projector wrapper (lazy loading singleton)
- Serialized models in `ml/models/`: `xgb_model.pkl`, `umap.pkl`

**Tests (`tests/`):**
- `conftest.py` - Pytest fixtures
- `test_health.py` - Health endpoint tests
- `test_predict.py` - Prediction endpoint tests
- `test_properties.py` - Properties endpoint tests

### Frontend (`frontend/`)
- **Framework:** React 19 + TypeScript + Vite
- **Entry point:** `src/main.tsx`

**Key Components (`src/components/`):**
- `MoleculeAnalyzer.tsx` - Main analysis component integrating all features
- `KetcherEditor.tsx` - Ketcher-based molecular structure editor
- `UmapPlot.tsx` - Plotly scatter plot for chemical space visualization
- `LipinskiProperties.tsx` - Lipinski Rule of 5 display
- `ToxicityExplanation.tsx` - Structural alerts and risk factors

**API Client (`src/api/client.ts`):**
- Functions: `apiPredict()`, `apiAnalyze()`, `apiGetUmapProjection()`, `apiGetLipinski()`, `apiExplainToxicity()`

## Data Flow

```
SMILES Input → Standardization → Mol Object → Feature Extraction (297 dims) → XGBoost Prediction
                                            └→ UMAP Projection (2D coordinates)
```

## Key Dependencies

- **Chemistry:** RDKit, scikit-fingerprints
- **ML:** XGBoost, scikit-learn, UMAP-learn
- **Frontend:** Ketcher (molecular editor), Plotly (visualization)
- **Testing:** pytest, httpx

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `CORS_ORIGINS` | Comma-separated allowed origins | `http://localhost:5173,http://localhost:3000,http://localhost:8000` |
| `VITE_API_BASE_URL` | Backend URL for frontend | - |

See `.env.example` for configuration template.
