"""
Simple test script to verify the project setup.

Run: python main.py
"""

import sys


def test_imports():
    """Test that all required packages can be imported."""
    print("Testing imports...")

    try:
        from rdkit import Chem
        print("  rdkit: OK")
    except ImportError as e:
        print(f"  rdkit: FAILED - {e}")
        return False

    try:
        import numpy as np
        print("  numpy: OK")
    except ImportError as e:
        print(f"  numpy: FAILED - {e}")
        return False

    try:
        import pandas as pd
        print("  pandas: OK")
    except ImportError as e:
        print(f"  pandas: FAILED - {e}")
        return False

    try:
        from skfp.fingerprints import MACCSFingerprint
        print("  scikit-fingerprints: OK")
    except ImportError as e:
        print(f"  scikit-fingerprints: FAILED - {e}")
        return False

    try:
        import umap
        print("  umap-learn: OK")
    except ImportError as e:
        print(f"  umap-learn: FAILED - {e}")
        return False

    try:
        from fastapi import FastAPI
        print("  fastapi: OK")
    except ImportError as e:
        print(f"  fastapi: FAILED - {e}")
        return False

    return True


def test_molecule_processing():
    """Test basic molecule processing."""
    print("\nTesting molecule processing...")

    from rdkit import Chem

    # Test SMILES parsing
    test_smiles = "CCO"  # ethanol
    mol = Chem.MolFromSmiles(test_smiles)

    if mol is None:
        print(f"  Failed to parse SMILES: {test_smiles}")
        return False

    print(f"  Parsed SMILES '{test_smiles}': OK")

    # Test MACCS fingerprint
    from skfp.fingerprints import MACCSFingerprint

    maccs_gen = MACCSFingerprint(n_jobs=1)
    fp = maccs_gen.transform([mol])
    print(f"  MACCS fingerprint shape: {fp.shape}")

    return True


def test_descriptors():
    """Test descriptor computation."""
    print("\nTesting descriptor computation...")

    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles("CCO")

    # Test some descriptors
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)

    print(f"  Molecular Weight: {mw:.2f}")
    print(f"  LogP: {logp:.2f}")

    return True


def main():
    print("=" * 50)
    print("Cardiotoxicity Prediction - Setup Test")
    print("=" * 50)

    all_passed = True

    if not test_imports():
        all_passed = False
        print("\nSome imports failed. Please install missing packages.")
    else:
        if not test_molecule_processing():
            all_passed = False
        if not test_descriptors():
            all_passed = False

    print("\n" + "=" * 50)
    if all_passed:
        print("All tests passed!")
        print("\nTo start the API server:")
        print("  uvicorn app:app --reload")
        print("\nTo train UMAP model:")
        print("  1. Place training_data.csv in project root")
        print("  2. Run umap_play.ipynb notebook")
    else:
        print("Some tests failed. Check the output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
