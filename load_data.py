from __future__ import annotations

import pandas as pd

from teiko.config import DATA_FILE, DB_FILE
from teiko.db import connect, initialize_database, load_dataframe
from teiko.validation import validate_input


def main() -> None:
    df = pd.read_csv(DATA_FILE, keep_default_na=False)
    report = validate_input(df)

    with connect(DB_FILE) as conn:
        initialize_database(conn)
        load_dataframe(conn, df)

    print(f"Loaded {report.rows} samples into {DB_FILE.name}")
    for warning in report.warnings:
        print(f"WARNING: {warning}")


if __name__ == "__main__":
    main()
