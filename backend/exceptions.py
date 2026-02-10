"""
Custom Exceptions

Standardized exceptions for the Cardiotoxicity Prediction API.
"""


class CardiotoxicityAPIError(Exception):
    """Base exception for all API errors."""

    def __init__(self, message: str, status_code: int = 400):
        self.message = message
        self.status_code = status_code
        super().__init__(self.message)


class InvalidSMILESError(CardiotoxicityAPIError):
    """Raised when SMILES string is invalid or cannot be parsed."""

    def __init__(self, message: str = "Invalid SMILES string"):
        super().__init__(message, status_code=400)


class MoleculeProcessingError(CardiotoxicityAPIError):
    """Raised when molecule processing fails."""

    def __init__(self, message: str = "Failed to process molecule"):
        super().__init__(message, status_code=400)


class ModelNotLoadedError(CardiotoxicityAPIError):
    """Raised when ML model is not available."""

    def __init__(self, message: str = "Model not loaded"):
        super().__init__(message, status_code=503)
