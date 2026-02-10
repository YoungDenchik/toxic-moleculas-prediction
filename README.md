# hERG Cardiotoxicity Predictor

A full-stack ML web application for predicting hERG cardiotoxicity of drug molecules. The hERG (human Ether-a-go-go-Related Gene) potassium channel is a critical target for drug safety assessment, as its inhibition can lead to potentially fatal cardiac arrhythmias.

## Features

- **Molecular Structure Input**: Draw molecules using the Ketcher molecular editor or input SMILES notation directly
- **Toxicity Prediction**: XGBoost classifier trained on ~12,000 compounds predicts toxic/non-toxic classification
- **Chemical Space Visualization**: Interactive UMAP plot showing where your molecule sits relative to known compounds
- **Lipinski Rule of 5**: Automatic calculation of drug-likeness properties
- **Structural Explanations**: Identification of structural alerts and risk factors associated with hERG toxicity

## Tech Stack

### Backend
- **Framework**: FastAPI + Uvicorn
- **Chemistry**: RDKit, scikit-fingerprints
- **ML**: XGBoost, scikit-learn, UMAP-learn
- **Testing**: pytest, httpx

### Frontend
- **Framework**: React 19 + TypeScript + Vite
- **Molecular Editor**: Ketcher
- **Visualization**: Plotly

## Quick Start

### Using Docker (Recommended)

**Development mode with hot reload:**
```bash
docker-compose -f docker-compose.dev.yml up
```
- Frontend: http://localhost:5173
- Backend: http://localhost:8000

**Production mode:**
```bash
docker-compose up
```
- Frontend: http://localhost:3000
- Backend: http://localhost:8000

### Local Development (without Docker)

**Backend:**
```bash
cd backend
pip install -r requirements.txt
uvicorn app:app --reload --host 0.0.0.0 --port 8000
```

**Frontend:**
```bash
cd frontend
npm install
npm run dev
```

### Running Tests

```bash
cd backend
pytest                        # Run all tests
pytest tests/test_predict.py  # Run specific test file
pytest -v                     # Verbose output
```

## Project Structure

```
├── backend/
│   ├── app.py                 # FastAPI entry point
│   ├── config.py              # Centralized configuration
│   ├── exceptions.py          # Custom API exceptions
│   ├── api/                   # API route handlers
│   │   ├── schemas.py         # Pydantic models
│   │   ├── predict.py         # /predict endpoint
│   │   ├── chemical_space.py  # /chemical-space/umap endpoint
│   │   ├── properties.py      # /properties/lipinski endpoint
│   │   ├── analyze.py         # /analyze (combined analysis)
│   │   └── explain.py         # /explain/toxicity endpoint
│   ├── services/              # Business logic layer
│   │   ├── molecule_service.py    # Centralized molecule processing
│   │   ├── prediction_pipeline.py # Prediction workflow
│   │   ├── umap_pipeline.py       # UMAP projection
│   │   ├── lipinski_pipeline.py   # Physicochemical properties
│   │   ├── explainability.py      # Structural alerts
│   │   ├── umap_data.py           # Training data loader
│   │   └── core_utils/            # Chemistry utilities
│   │       ├── standardization.py # SMILES canonicalization
│   │       ├── mol.py             # RDKit molecule creation
│   │       ├── features.py        # Feature extraction
│   │       └── descriptors.py     # Descriptor configuration
│   ├── ml/                    # ML models
│   │   ├── model.py           # XGBoost wrapper
│   │   ├── umap_model.py      # UMAP wrapper
│   │   └── models/            # Serialized model files
│   └── tests/                 # Test suite
│       ├── conftest.py        # Pytest fixtures
│       ├── test_health.py     # Health endpoint tests
│       ├── test_predict.py    # Prediction tests
│       └── test_properties.py # Properties tests
├── frontend/
│   ├── src/
│   │   ├── main.tsx           # React entry point
│   │   ├── types.ts           # TypeScript interfaces
│   │   ├── components/        # UI components
│   │   │   ├── MoleculeAnalyzer.tsx
│   │   │   ├── KetcherEditor.tsx
│   │   │   ├── UmapPlot.tsx
│   │   │   ├── LipinskiProperties.tsx
│   │   │   └── ToxicityExplanation.tsx
│   │   └── api/
│   │       └── client.ts      # API client functions
│   └── package.json
├── docker-compose.yml         # Production Docker config
├── docker-compose.dev.yml     # Development Docker config
└── .env.example               # Environment variables template
```

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Health check |
| `/health` | GET | Health status |
| `/predict` | POST | Predict toxicity from SMILES |
| `/chemical-space/umap` | POST | Get 2D UMAP coordinates |
| `/properties/lipinski` | POST | Calculate Lipinski properties |
| `/analyze` | POST | Full analysis (all above combined) |
| `/explain/toxicity` | POST | Get structural alerts and risk factors |
| `/umap` | GET | Get training data UMAP coordinates |

## ML Pipeline

The prediction pipeline follows these steps:

1. **SMILES Standardization**: Canonicalize input molecular structure
2. **Feature Extraction**: Generate 297-dimensional feature vector
   - 167 MACCS fingerprint bits
   - 130 RDKit molecular descriptors
3. **Prediction**: XGBoost classifier outputs probability and class
4. **Visualization**: UMAP projects molecule into 2D chemical space

## Dataset

The model was trained on a curated dataset of ~12,000 compounds with experimentally determined hERG activity. The dataset is approximately balanced (53% toxic, 47% non-toxic) and was split using scaffold-based splitting to ensure structural diversity between train and validation sets.

## Configuration

Copy `.env.example` to `.env` and configure as needed:

| Variable | Description | Default |
|----------|-------------|---------|
| `CORS_ORIGINS` | Comma-separated allowed origins | `http://localhost:5173,http://localhost:3000,http://localhost:8000` |
| `VITE_API_BASE_URL` | Backend API URL for frontend | - |

## License

MIT
