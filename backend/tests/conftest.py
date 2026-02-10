"""
Pytest configuration and fixtures.
"""

import pytest
from fastapi.testclient import TestClient

from backend.app import app


@pytest.fixture
def client():
    """Create a test client for the FastAPI app."""
    return TestClient(app)


@pytest.fixture
def valid_smiles():
    """Return a valid SMILES string for testing."""
    return "CCO"  # Ethanol


@pytest.fixture
def invalid_smiles():
    """Return an invalid SMILES string for testing."""
    return "not_a_valid_smiles_string"


@pytest.fixture
def toxic_smiles():
    """Return a SMILES string likely to be predicted as toxic."""
    # Dofetilide - known hERG blocker
    return "CC(C)N1CCN(CCOc2ccc(NS(=O)(=O)c3ccc(C)cc3)cc2)CC1"


@pytest.fixture
def non_toxic_smiles():
    """Return a simple molecule likely to be predicted as non-toxic."""
    return "C"  # Methane
