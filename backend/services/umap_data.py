"""
UMAP Training Data Service

Loads and caches precomputed UMAP coordinates for training data.
"""

import csv
from functools import lru_cache
from pathlib import Path
from typing import List, Dict, Any


@lru_cache(maxsize=1)
def load_umap_training_data() -> List[Dict[str, Any]]:
    """
    Load precomputed UMAP coordinates from CSV.

    Returns list of dicts matching frontend UmapPoint structure:
    - id: training_{index}
    - x: UMAP_1 coordinate
    - y: UMAP_2 coordinate
    - label: 1 if Toxic, 0 if Non-toxic
    - smiles: molecule SMILES string
    """
    csv_path = Path(__file__).parent.parent / "ml" / "models" / "umap_coords.csv"

    points = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader):
            points.append({
                "id": f"training_{idx}",
                "x": float(row["UMAP_1"]),
                "y": float(row["UMAP_2"]),
                "label": 1 if row["Activity"] == "Toxic" else 0,
                "smiles": row["SMILES"],
            })

    return points
