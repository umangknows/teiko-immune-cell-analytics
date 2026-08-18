from __future__ import annotations

import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import CELL_POPULATIONS


def exploratory_prediction(freq: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide = (
        freq.pivot_table(
            index=["subject", "sample", "response"],
            columns="population",
            values="percentage",
            aggfunc="first",
        )
        .reset_index()
        .dropna(subset=CELL_POPULATIONS)
    )
    wide["target"] = (wide["response"] == "yes").astype(int)
    X = wide[CELL_POPULATIONS]
    y = wide["target"]

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    model = Pipeline(
        [
            ("scale", StandardScaler()),
            ("logit", LogisticRegression(max_iter=1000, random_state=42)),
        ]
    )
    probabilities = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
    auc = roc_auc_score(y, probabilities)

    model.fit(X, y)
    coefficients = model.named_steps["logit"].coef_[0]
    importance = pd.DataFrame(
        {
            "population": CELL_POPULATIONS,
            "logistic_coefficient": coefficients,
            "absolute_importance": abs(coefficients),
        }
    ).sort_values("absolute_importance", ascending=False)
    summary = pd.DataFrame(
        {
            "model": ["logistic_regression_baseline_pbmc"],
            "samples": [len(wide)],
            "cv_folds": [5],
            "cross_validated_auc": [round(float(auc), 4)],
            "note": ["Exploratory signal only; not clinically validated."],
        }
    )
    return summary, importance

