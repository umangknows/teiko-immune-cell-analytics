import pandas as pd
import pytest

from teiko.validation import validate_input


def valid_row(**overrides):
    row = {
        "project": "prj1",
        "subject": "sbj1",
        "condition": "melanoma",
        "age": 50,
        "sex": "M",
        "treatment": "miraclib",
        "response": "yes",
        "sample": "sample1",
        "sample_type": "PBMC",
        "time_from_treatment_start": 0,
        "b_cell": 10,
        "cd8_t_cell": 20,
        "cd4_t_cell": 30,
        "nk_cell": 15,
        "monocyte": 25,
    }
    row.update(overrides)
    return row


def test_validation_accepts_valid_rows():
    report = validate_input(pd.DataFrame([valid_row()]))
    assert report.rows == 1


def test_validation_rejects_duplicate_samples():
    df = pd.DataFrame([valid_row(), valid_row(subject="sbj2")])
    with pytest.raises(ValueError, match="Duplicate sample"):
        validate_input(df)


def test_validation_rejects_negative_counts():
    df = pd.DataFrame([valid_row(b_cell=-1)])
    with pytest.raises(ValueError, match="negative"):
        validate_input(df)

