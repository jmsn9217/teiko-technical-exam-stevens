"""Initialize the Teiko SQLite database and load cell-count.csv.

Run from the repo root: ``python load_data.py``
"""

from pathlib import Path
import sqlite3
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "teiko.db"
CSV_PATH = ROOT / "cell-count.csv"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def create_schema(conn: sqlite3.Connection) -> None:
    cur = conn.cursor()
    cur.executescript(
        """
        DROP TABLE IF EXISTS cell_counts;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS subjects;
        DROP TABLE IF EXISTS projects;

        CREATE TABLE projects (
            project_id TEXT PRIMARY KEY
        );

        CREATE TABLE subjects (
            subject_id TEXT PRIMARY KEY,
            project_id TEXT NOT NULL REFERENCES projects(project_id),
            condition  TEXT NOT NULL,
            age        INTEGER,
            sex        TEXT,
            treatment  TEXT,
            response   TEXT
        );

        CREATE TABLE samples (
            sample_id                  TEXT PRIMARY KEY,
            subject_id                 TEXT NOT NULL REFERENCES subjects(subject_id),
            sample_type                TEXT NOT NULL,
            time_from_treatment_start  INTEGER
        );

        CREATE TABLE cell_counts (
            sample_id  TEXT NOT NULL REFERENCES samples(sample_id),
            population TEXT NOT NULL,
            count      INTEGER NOT NULL,
            PRIMARY KEY (sample_id, population)
        );

        CREATE INDEX idx_samples_subject       ON samples(subject_id);
        CREATE INDEX idx_subjects_project      ON subjects(project_id);
        CREATE INDEX idx_subjects_cond_treat   ON subjects(condition, treatment);
        CREATE INDEX idx_cell_counts_pop       ON cell_counts(population);
        """
    )


def load_csv(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = df.rename(columns={"project": "project_id", "subject": "subject_id", "sample": "sample_id"})
    return df


def insert_projects(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    rows = [(pid,) for pid in sorted(df["project_id"].unique())]
    conn.executemany("INSERT INTO projects (project_id) VALUES (?)", rows)


def insert_subjects(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    cols = ["subject_id", "project_id", "condition", "age", "sex", "treatment", "response"]
    subjects = df[cols].drop_duplicates(subset="subject_id").where(pd.notna(df), None)
    # pandas .where(notna, None) keeps NaN floats; convert NaN -> None explicitly for sqlite
    rows = [tuple(None if pd.isna(v) else v for v in r) for r in subjects.itertuples(index=False, name=None)]
    conn.executemany(
        "INSERT INTO subjects (subject_id, project_id, condition, age, sex, treatment, response) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        rows,
    )


def insert_samples(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    cols = ["sample_id", "subject_id", "sample_type", "time_from_treatment_start"]
    samples = df[cols].drop_duplicates(subset="sample_id")
    rows = [tuple(None if pd.isna(v) else v for v in r) for r in samples.itertuples(index=False, name=None)]
    conn.executemany(
        "INSERT INTO samples (sample_id, subject_id, sample_type, time_from_treatment_start) "
        "VALUES (?, ?, ?, ?)",
        rows,
    )


def insert_cell_counts(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    long = df.melt(
        id_vars=["sample_id"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    rows = list(long.itertuples(index=False, name=None))
    conn.executemany(
        "INSERT INTO cell_counts (sample_id, population, count) VALUES (?, ?, ?)",
        rows,
    )


def main() -> int:
    if not CSV_PATH.exists():
        print(f"ERROR: {CSV_PATH} not found", file=sys.stderr)
        return 1

    print(f"Loading {CSV_PATH.name} -> {DB_PATH.name}")
    df = load_csv(CSV_PATH)

    if DB_PATH.exists():
        DB_PATH.unlink()

    with sqlite3.connect(DB_PATH) as conn:
        conn.execute("PRAGMA foreign_keys = ON")
        create_schema(conn)
        insert_projects(conn, df)
        insert_subjects(conn, df)
        insert_samples(conn, df)
        insert_cell_counts(conn, df)
        conn.commit()

        for table in ("projects", "subjects", "samples", "cell_counts"):
            (n,) = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()
            print(f"  {table:<12} {n:>6} rows")

    print(f"Done. Database written to {DB_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
