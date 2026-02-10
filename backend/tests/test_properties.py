"""
Tests for properties endpoints.
"""


class TestLipinskiEndpoint:
    """Tests for the /properties/lipinski endpoint."""

    def test_lipinski_valid_smiles(self, client, valid_smiles):
        """Test Lipinski properties with a valid SMILES string."""
        response = client.post(
            "/properties/lipinski", json={"smiles": valid_smiles}
        )
        assert response.status_code == 200
        data = response.json()
        assert "mw" in data
        assert "logp" in data
        assert "hbd" in data
        assert "hba" in data
        assert "ro5_violations" in data

    def test_lipinski_invalid_smiles(self, client, invalid_smiles):
        """Test Lipinski properties with an invalid SMILES string."""
        response = client.post(
            "/properties/lipinski", json={"smiles": invalid_smiles}
        )
        assert response.status_code == 400


class TestExplainEndpoint:
    """Tests for the /explain/toxicity endpoint."""

    def test_explain_valid_smiles(self, client, valid_smiles):
        """Test toxicity explanation with a valid SMILES string."""
        response = client.post(
            "/explain/toxicity", json={"smiles": valid_smiles}
        )
        assert response.status_code == 200
        data = response.json()
        assert "structural_alerts" in data
        assert "risk_factors" in data
        assert isinstance(data["structural_alerts"], list)
        assert isinstance(data["risk_factors"], list)

    def test_explain_invalid_smiles(self, client, invalid_smiles):
        """Test toxicity explanation with an invalid SMILES string."""
        response = client.post(
            "/explain/toxicity", json={"smiles": invalid_smiles}
        )
        assert response.status_code == 400
