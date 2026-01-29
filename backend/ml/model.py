import joblib
import numpy as np



class CardiotoxicityModel:
    def __init__(self, model_path: str):
        self.package = joblib.load(model_path)
        self.model = self.package["model"]

    def predict_proba(self, X: np.ndarray) -> float:
        """
        Expects shape (n_features,)
        """
        X = X.reshape(1, -1)
        return float(self.model.predict_proba(X)[0, 1])


# singleton
model = CardiotoxicityModel("backend/ml/models/xgb_model.pkl")
