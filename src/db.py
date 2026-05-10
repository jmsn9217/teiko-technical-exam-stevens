"""Shared SQL queries used by analysis.py and the Streamlit dashboard.

Keeping the queries here (instead of inline) means analysis.py and the dashboard
read the data the same way, so a fix to a query is a fix everywhere.
"""

from pathlib import Path
import sqlite3

import pandas as pd

DB_PATH = Path(__file__).resolve().parent.parent / "teiko.db"


def connect(db_path: Path | str = DB_PATH) -> sqlite3.Connection:
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def summary_table(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 2: long-format relative frequency per (sample, population).

    Columns: sample, total_count, population, count, percentage.
    Percentage is computed in SQL with a window function so the math lives
    next to the data.
    """
    sql = """
        SELECT
            cc.sample_id                                                 AS sample,
            SUM(cc.count) OVER (PARTITION BY cc.sample_id)               AS total_count,
            cc.population                                                AS population,
            cc.count                                                     AS count,
            100.0 * cc.count
                / SUM(cc.count) OVER (PARTITION BY cc.sample_id)         AS percentage
        FROM cell_counts cc
        ORDER BY cc.sample_id, cc.population
    """
    return pd.read_sql_query(sql, conn)


def responder_comparison_frame(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 3 input: per-sample relative frequencies for melanoma + miraclib + PBMC,
    annotated with response. Long format (one row per sample-population).
    """
    sql = """
        WITH per_sample AS (
            SELECT
                cc.sample_id,
                cc.population,
                cc.count,
                SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count
            FROM cell_counts cc
        )
        SELECT
            ps.sample_id                              AS sample,
            ps.population                             AS population,
            ps.count                                  AS count,
            ps.total_count                            AS total_count,
            100.0 * ps.count / ps.total_count         AS percentage,
            sub.response                              AS response,
            sub.subject_id                            AS subject_id
        FROM per_sample ps
        JOIN samples  s   ON s.sample_id  = ps.sample_id
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition  = 'melanoma'
          AND sub.treatment  = 'miraclib'
          AND s.sample_type  = 'PBMC'
          AND sub.response  IN ('yes', 'no')
        ORDER BY ps.sample_id, ps.population
    """
    return pd.read_sql_query(sql, conn)


def baseline_miraclib_melanoma_pbmc_samples(conn: sqlite3.Connection) -> pd.DataFrame:
    """Part 4: melanoma PBMC samples at time=0 from miraclib-treated subjects."""
    sql = """
        SELECT
            s.sample_id,
            s.subject_id,
            sub.project_id,
            sub.response,
            sub.sex,
            s.time_from_treatment_start
        FROM samples s
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
          AND s.time_from_treatment_start = 0
        ORDER BY s.sample_id
    """
    return pd.read_sql_query(sql, conn)


def b_cell_count_for_samples(conn: sqlite3.Connection, sample_ids: list[str]) -> pd.DataFrame:
    """B-cell counts for a given list of sample ids (Part 4 follow-up)."""
    if not sample_ids:
        return pd.DataFrame(columns=["sample_id", "count"])
    placeholders = ",".join("?" for _ in sample_ids)
    sql = f"""
        SELECT sample_id, count
        FROM cell_counts
        WHERE population = 'b_cell'
          AND sample_id IN ({placeholders})
    """
    return pd.read_sql_query(sql, conn, params=sample_ids)
