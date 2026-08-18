from __future__ import annotations

import pandas as pd

from .config import DB_FILE
from .db import read_sql


def frequency_summary(db_path=DB_FILE) -> pd.DataFrame:
    return read_sql(
        """
        SELECT sample, total_count, population, count, percentage
        FROM sample_cell_frequencies
        ORDER BY sample, population
        """,
        db_path,
    )


def frequency_summary_with_metadata(db_path=DB_FILE) -> pd.DataFrame:
    return read_sql(
        """
        SELECT
            project,
            subject,
            condition,
            age,
            sex,
            treatment,
            response,
            sample,
            sample_type,
            time_from_treatment_start,
            total_count,
            population,
            count,
            percentage
        FROM sample_cell_frequencies
        ORDER BY sample, population
        """,
        db_path,
    )


def sample_metadata(db_path=DB_FILE) -> pd.DataFrame:
    return read_sql(
        """
        SELECT project, subject, condition, age, sex, treatment, response,
               sample, sample_type, time_from_treatment_start
        FROM sample_metadata
        ORDER BY sample
        """,
        db_path,
    )


def responder_frequency_data(db_path=DB_FILE, baseline_only: bool = True) -> pd.DataFrame:
    time_filter = "AND time_from_treatment_start = 0" if baseline_only else ""
    return read_sql(
        f"""
        SELECT
            project,
            subject,
            sample,
            time_from_treatment_start,
            response,
            population,
            percentage
        FROM sample_cell_frequencies
        WHERE condition = 'melanoma'
          AND treatment = 'miraclib'
          AND sample_type = 'PBMC'
          AND response IN ('yes', 'no')
          {time_filter}
        ORDER BY population, response, subject, time_from_treatment_start
        """,
        db_path,
    )


def baseline_subset_summary(db_path=DB_FILE) -> dict[str, pd.DataFrame]:
    samples = read_sql(
        """
        SELECT DISTINCT project, subject, response, sex, sample
        FROM sample_metadata
        WHERE condition = 'melanoma'
          AND treatment = 'miraclib'
          AND sample_type = 'PBMC'
          AND time_from_treatment_start = 0
        """,
        db_path,
    )
    return {
        "samples_by_project": samples.groupby("project", as_index=False)["sample"]
        .nunique()
        .rename(columns={"sample": "sample_count"}),
        "subjects_by_response": samples.groupby("response", as_index=False)["subject"]
        .nunique()
        .rename(columns={"subject": "subject_count"}),
        "subjects_by_sex": samples.groupby("sex", as_index=False)["subject"]
        .nunique()
        .rename(columns={"subject": "subject_count"}),
    }


def melanoma_male_responder_bcell_average(db_path=DB_FILE) -> float:
    result = read_sql(
        """
        SELECT AVG(cc.count) AS avg_b_cell
        FROM cell_counts cc
        JOIN cell_populations cp ON cp.population_id = cc.population_id
        JOIN sample_metadata sm ON sm.sample = cc.sample_id
        WHERE cp.population_name = 'b_cell'
          AND sm.condition = 'melanoma'
          AND sm.sex = 'M'
          AND sm.response = 'yes'
          AND sm.time_from_treatment_start = 0
        """,
        db_path,
    )
    return round(float(result.loc[0, "avg_b_cell"]), 2)


def dataset_overview(db_path=DB_FILE) -> dict[str, int]:
    overview = read_sql(
        """
        SELECT
            (SELECT COUNT(*) FROM projects) AS projects,
            (SELECT COUNT(*) FROM subjects) AS subjects,
            (SELECT COUNT(*) FROM samples) AS samples,
            (SELECT COUNT(*) FROM cell_counts) AS cell_count_records
        """,
        db_path,
    )
    return {key: int(value) for key, value in overview.iloc[0].to_dict().items()}
