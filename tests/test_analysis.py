from pathlib import Path

import pandas as pd

from teiko.config import OUTPUT_DIR
from teiko.analysis import frequency_summary, melanoma_male_responder_bcell_average
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
            }
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
    assert round(summary["percentage"].sum(), 6) == 100


def test_average_bcell_query_formats_expected_value():
    db_path = build_test_db()
    assert melanoma_male_responder_bcell_average(db_path) == 10.00
