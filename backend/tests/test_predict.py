"""
Tests for prediction endpoints.
"""

import pytest


def test_predict_valid_smiles(client, valid_smiles):
    """Test prediction with a valid SMILES string."""
    response = client.post("/predict", json={"smiles": valid_smiles})
    assert response.status_code == 200
    data = response.json()
    assert "canonical_smiles" in data
    assert "prediction" in data
    assert data["prediction"]["label"] in ["toxic", "non-toxic"]
    assert 0 <= data["prediction"]["probability"] <= 1


def test_predict_invalid_smiles(client, invalid_smiles):
    """Test prediction with an invalid SMILES string."""
    response = client.post("/predict", json={"smiles": invalid_smiles})
    assert response.status_code == 400
    data = response.json()
    assert "detail" in data


def test_predict_empty_smiles(client):
    """Test prediction with an empty SMILES string."""
    response = client.post("/predict", json={"smiles": ""})
    assert response.status_code == 400


def test_predict_missing_smiles(client):
    """Test prediction with missing SMILES field."""
    response = client.post("/predict", json={})
    assert response.status_code == 422  # Validation error


class TestAnalyzeEndpoint:
    """Tests for the /analyze endpoint."""

    def test_analyze_valid_smiles(self, client, valid_smiles):
        """Test full analysis with a valid SMILES string."""
        response = client.post("/analyze", json={"smiles": valid_smiles})
        assert response.status_code == 200
        data = response.json()
        assert "canonical_smiles" in data
        assert "prediction" in data
        assert "umap" in data
        assert "lipinski" in data

    def test_analyze_invalid_smiles(self, client, invalid_smiles):
        """Test analysis with an invalid SMILES string."""
        response = client.post("/analyze", json={"smiles": invalid_smiles})
        assert response.status_code == 400


class TestChemicalSpaceEndpoint:
    """Tests for the /chemical-space/umap endpoint."""

    def test_umap_projection_valid_smiles(self, client, valid_smiles):
        """Test UMAP projection with a valid SMILES string."""
        response = client.post(
            "/chemical-space/umap", json={"smiles": valid_smiles}
        )
        assert response.status_code == 200
        data = response.json()
        assert "x" in data
        assert "y" in data
        assert isinstance(data["x"], (int, float))
        assert isinstance(data["y"], (int, float))

    def test_umap_training_data(self, client):
        """Test fetching UMAP training data."""
        response = client.get("/umap")
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        if len(data) > 0:
            point = data[0]
            assert "id" in point
            assert "x" in point
            assert "y" in point
            assert "label" in point
