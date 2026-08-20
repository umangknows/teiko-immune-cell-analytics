from pathlib import Path

import pandas as pd

from teiko.config import OUTPUT_DIR
from teiko.analysis import (
    baseline_subset_summary,
    frequency_summary,
    melanoma_male_responder_bcell_average,
)
from teiko.db import connect, initialize_database, load_dataframe
from teiko.validation import validate_input


def build_test_db() -> Path:
    df = pd.DataFrame(
        [
            {
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
            },
            {
                "project": "prj1",
                "subject": "sbj2",
                "condition": "melanoma",
                "age": 61,
                "sex": "F",
                "treatment": "miraclib",
                "response": "no",
                "sample": "sample2",
                "sample_type": "PBMC",
                "time_from_treatment_start": 0,
                "b_cell": 12,
                "cd8_t_cell": 18,
                "cd4_t_cell": 35,
                "nk_cell": 16,
                "monocyte": 19,
            },
            {
                "project": "prj2",
                "subject": "sbj3",
                "condition": "melanoma",
                "age": 58,
                "sex": "M",
                "treatment": "phauximab",
                "response": "yes",
                "sample": "sample3",
                "sample_type": "WB",
                "time_from_treatment_start": 0,
                "b_cell": 30,
                "cd8_t_cell": 10,
                "cd4_t_cell": 20,
                "nk_cell": 15,
                "monocyte": 25,
            },
            {
                "project": "prj2",
                "subject": "sbj4",
                "condition": "melanoma",
                "age": 42,
                "sex": "M",
                "treatment": "miraclib",
                "response": "yes",
                "sample": "sample4",
                "sample_type": "WB",
                "time_from_treatment_start": 0,
                "b_cell": 50,
                "cd8_t_cell": 10,
                "cd4_t_cell": 20,
                "nk_cell": 10,
                "monocyte": 10,
            },
        ]
    )
    validate_input(df)
    OUTPUT_DIR.mkdir(exist_ok=True)
    db_path = OUTPUT_DIR / "test_analysis.db"
    with connect(db_path) as conn:
        initialize_database(conn)
        load_dataframe(conn, df)
    return db_path


def test_frequency_percentages_sum_to_100():
    db_path = build_test_db()
    summary = frequency_summary(db_path)
    assert list(summary.columns) == ["sample", "total_count", "population", "count", "percentage"]
    assert summary.shape[0] == 20
    assert summary.groupby("sample")["percentage"].sum().round(6).eq(100).all()


def test_average_bcell_query_uses_all_sample_and_treatment_types():
    db_path = build_test_db()
    assert melanoma_male_responder_bcell_average(db_path) == 30.00


def test_baseline_subset_summary_uses_miraclib_pbmc_baseline_only():
    db_path = build_test_db()
    subset = baseline_subset_summary(db_path)

    assert subset["baseline_melanoma_pbmc_miraclib_samples"].to_dict("records") == [
        {
            "project": "prj1",
            "subject": "sbj1",
            "response": "yes",
            "sex": "M",
            "sample": "sample1",
        },
        {
            "project": "prj1",
            "subject": "sbj2",
            "response": "no",
            "sex": "F",
            "sample": "sample2",
        },
    ]
    assert subset["samples_by_project"].to_dict("records") == [
        {"project": "prj1", "sample_count": 2}
    ]
    assert subset["subjects_by_response"].to_dict("records") == [
        {"response": "no", "subject_count": 1},
        {"response": "yes", "subject_count": 1},
    ]
    assert subset["subjects_by_sex"].to_dict("records") == [
        {"sex": "F", "subject_count": 1},
        {"sex": "M", "subject_count": 1},
    ]
