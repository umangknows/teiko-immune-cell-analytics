from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd

from .config import CELL_POPULATIONS, DB_FILE
from .schema import SCHEMA_SQL


def connect(db_path: Path = DB_FILE) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def initialize_database(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA_SQL)


def load_dataframe(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    projects = [(project,) for project in sorted(df["project"].unique())]
    conn.executemany("INSERT INTO projects(project_id) VALUES (?)", projects)

    subjects = (
        df[["subject", "project", "condition", "age", "sex", "treatment", "response"]]
        .drop_duplicates("subject")
        .assign(response=lambda x: x["response"].replace("", None))
    )
    conn.executemany(
        """
        INSERT INTO subjects(subject_id, project_id, condition, age, sex, treatment, response)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        subjects.itertuples(index=False, name=None),
    )

    samples = df[["sample", "subject", "sample_type", "time_from_treatment_start"]]
    conn.executemany(
        """
        INSERT INTO samples(sample_id, subject_id, sample_type, time_from_treatment_start)
        VALUES (?, ?, ?, ?)
        """,
        samples.itertuples(index=False, name=None),
    )

    populations = [(i + 1, name) for i, name in enumerate(CELL_POPULATIONS)]
    conn.executemany(
        "INSERT INTO cell_populations(population_id, population_name) VALUES (?, ?)",
        populations,
    )
    population_ids = {name: i + 1 for i, name in enumerate(CELL_POPULATIONS)}

    count_rows = []
    for row in df.itertuples(index=False):
        sample = getattr(row, "sample")
        for population in CELL_POPULATIONS:
            count_rows.append((sample, population_ids[population], int(getattr(row, population))))

    conn.executemany(
        "INSERT INTO cell_counts(sample_id, population_id, count) VALUES (?, ?, ?)",
        count_rows,
    )
    conn.commit()


def read_sql(query: str, db_path: Path = DB_FILE, params: tuple = ()) -> pd.DataFrame:
    with connect(db_path) as conn:
        return pd.read_sql_query(query, conn, params=params)

