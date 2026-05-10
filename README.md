# Teiko technical exam — Josh Stevens

Analysis pipeline + interactive dashboard for cell-count data from Bob Loblaw's
clinical trial of miraclib. Built around a small SQLite database, a pandas /
scipy analysis layer, and a Streamlit dashboard.

## Quick start (Codespaces or local)

```bash
make setup       # installs requirements.txt
make pipeline    # builds teiko.db, generates everything in outputs/
make dashboard   # serves the Streamlit app on http://localhost:8501
```

`make pipeline` runs `load_data.py` (Part 1) and then `python -m src.analysis`
(Parts 2-4) end-to-end, with no manual intervention.

## Dashboard

**Live link:** _to be deployed to Streamlit Community Cloud after submission._
Run locally with `make dashboard`.

The dashboard has three tabs, one per analytical part:

1. **Per-sample frequencies** — interactive view + download of the Part 2 summary table.
2. **Responders vs non-responders** — Plotly boxplot + the Mann-Whitney U / Bonferroni stats.
3. **Baseline subset** — counts and the average B-cell metric from Part 4.

## Repository layout

```
load_data.py          Part 1: schema + CSV load -> teiko.db (root, runnable directly)
cell-count.csv        input data (committed for reproducibility)
src/
  db.py               shared SQL queries used by analysis + dashboard
  analysis.py         Parts 2-4 orchestrator (writes to outputs/)
dashboard/
  app.py              Streamlit dashboard (reads teiko.db via src/db.py)
outputs/              generated tables (CSV) and plots (PNG, JSON)
Makefile              setup / pipeline / dashboard targets
requirements.txt
```

### Design choices

- **One source of truth for queries.** Both `analysis.py` and the dashboard
  call into `src/db.py`. If a query is wrong, fixing it in one place fixes it
  everywhere — no chance of drift between the static report and the dashboard.
- **Computed columns in SQL, not Python.** `total_count` and `percentage` are
  computed at query time with a window function, so they're always consistent
  with whatever `count` rows exist. No risk of stale derived data.
- **Long-format `cell_counts`.** Adding a sixth cell population is a row, not a
  schema change — see schema rationale below.
- **Headless matplotlib.** `matplotlib.use("Agg")` so the pipeline runs cleanly
  in CI / Codespaces with no display.

## Database schema

```
projects(project_id PK)
subjects(subject_id PK, project_id FK, condition, age, sex, treatment, response)
samples(sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
cell_counts(sample_id FK, population, count, PRIMARY KEY (sample_id, population))
```

## Database Design Rationale
Note: the prompt references 'indication' but the CSV column is 'condition'. This schema uses 'condition' to match the source data.

Cell counts get their own table with the PK being (sample_id, population), because this (long format) will make it easier to scale.
Having one row per population sample means that adding or removing cell populations is as easy as adding or removing a row rather than a column, i.e. wouldn't require a schema change. 
It also allows for easier analysis, e.g. querying with GROUP BY.
The total_count and percentage should be computed at query time to ensure that the data is fresh.

Separate domains for samples and subjects. 
Ensures easy and accurate analysis of subject specific data. 
If subject attributes were placed in each sample row, it would make things more complicated, like querying for responders vs. non-responders. 
This way is more efficient as it eliminates the need for things like DISTINCT in queries.
Normalized to 3NF. 
Subject attributes don't vary within a subject, so they belong on the subject entity.
This prevents inconsistency at scale.
However, if you wanted to centralize subject data across projects, you could add a subject_attributes table (attributes_id PK, subject FK) with subject attributes that could vary depending on the project, but still map them back to the original subject record. Certain static attributes like date_of_birth, sex, etc could then be extracted and stored only on the subject record.

Projects are a separate entity as well.
Having projects as its own entity makes it easy to add project metadata (sponsor, start date, principal investigator).
It also allows for the potential to merge and thus scale data more easily, for example to have one centralized table for all projects where multiple subjects could be mapped to multiple projects.

Indexes to add at scale: 
- samples(subject_id): supports joining sample data to subject attributes for any subject-level filtering (Parts 3 and 4 both need this).
- subjects(project_id, condition, treatment): supports the Part 4 cohort filter (melanoma + miraclib + project breakdown).
- cell_counts(population): supports per-population aggregation across samples (Parts 2 and 3).

This schema scales linearly. 
The queries from Parts 2-4 become easy joins. 
Aggregations, like relative frequency by sample, can be done via functions or views, which is great when the data scales up.

## Code structure overview

| File | Responsibility |
|---|---|
| `load_data.py` | Build the schema, load the CSV. Idempotent (drops + recreates). |
| `src/db.py` | All SQL queries — one place to change a query. |
| `src/analysis.py` | Part 2 (summary), Part 3 (stats + boxplot), Part 4 (subset). Writes to `outputs/`. |
| `dashboard/app.py` | Streamlit UI. Read-only — never writes to the DB. |

The split is deliberate: writes go through `load_data.py` only, the analysis
script and the dashboard are pure readers, and they share their query layer.
That makes each piece testable in isolation and the dashboard safe to deploy
publicly.

## Statistical method (Part 3)

Mann-Whitney U test (`scipy.stats.mannwhitneyu`, two-sided) per population,
with a Bonferroni correction across the 5 populations.

- **Why Mann-Whitney, not t-test:** cell-frequency distributions are typically
  right-skewed and small-sample, so the t-test's normality assumption is
  shaky. Mann-Whitney makes no distributional assumption.
- **Why Bonferroni:** five simultaneous tests inflate the family-wise error
  rate; Bonferroni is the most conservative correction and the cleanest one
  to defend ("if it's significant after Bonferroni, it's really significant").
  For a more powerful but slightly less conservative test, swap to
  Benjamini-Hochberg FDR — one-line change.

In the current data, no population reaches Bonferroni-corrected p < 0.05;
`cd4_t_cell` (raw p ≈ 0.013) and `b_cell` (raw p ≈ 0.056) are the most
suggestive but don't survive correction. That's an honest finding worth
reporting to Yah, not a reason to fish for significance.

## Outputs (generated by `make pipeline`)

```
outputs/
  part2_summary_table.csv          per-sample relative frequencies, long format
  part3_responder_frequencies.csv  Part 3 input frame (melanoma+miraclib+PBMC)
  part3_stats.csv                  Mann-Whitney U + Bonferroni results
  part3_boxplot.png                static boxplot for the report
  part4_baseline_samples.csv       baseline cohort listing
  part4_summary.json               counts + average B-cell metric
```

## Reproducing in Codespaces

```bash
make setup
make pipeline       # generates teiko.db and everything in outputs/
make dashboard      # opens http://localhost:8501
```
