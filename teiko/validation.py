from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from .config import CELL_POPULATIONS, REQUIRED_COLUMNS


@dataclass(frozen=True)
class ValidationReport:
    rows: int
    warnings: list[str]


def validate_input(df: pd.DataFrame) -> ValidationReport:
    missing_columns = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    warnings: list[str] = []
    work = df[REQUIRED_COLUMNS].copy()

    if work["sample"].duplicated().any():
        duplicates = work.loc[work["sample"].duplicated(), "sample"].head(5).tolist()
        raise ValueError(f"Duplicate sample identifiers found: {duplicates}")

    subject_time_dupes = work.duplicated(
        subset=["subject", "sample_type", "time_from_treatment_start"]
    )
    if subject_time_dupes.any():
        examples = work.loc[
            subject_time_dupes, ["subject", "sample_type", "time_from_treatment_start"]
        ].head(5)
        raise ValueError(f"Duplicate subject/sample_type/time rows found: {examples.to_dict('records')}")

    for col in ["age", "time_from_treatment_start", *CELL_POPULATIONS]:
        numeric = pd.to_numeric(work[col], errors="coerce")
        if numeric.isna().any():
            raise ValueError(f"Column {col} contains non-numeric values")
        if col in CELL_POPULATIONS and (numeric < 0).any():
            raise ValueError(f"Column {col} contains negative cell counts")

    totals = work[CELL_POPULATIONS].apply(pd.to_numeric, errors="raise").sum(axis=1)
    if (totals <= 0).any():
        raise ValueError("At least one sample has a non-positive total cell count")

    expected = {
        "sex": {"M", "F"},
        "sample_type": {"PBMC", "WB"},
        "response": {"yes", "no", ""},
    }
    for col, allowed in expected.items():
        observed = set(work[col].fillna("").astype(str))
        unexpected = sorted(observed - allowed)
        if unexpected:
            warnings.append(f"{col} contains unexpected values: {unexpected}")

    return ValidationReport(rows=len(work), warnings=warnings)

