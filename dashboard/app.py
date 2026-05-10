"""Streamlit dashboard for the Teiko technical exam.

Run with:  streamlit run dashboard/app.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import plotly.express as px  # type: ignore
import streamlit as st

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src import db  # noqa: E402

POPULATIONS = ["b_cell", "cd4_t_cell", "cd8_t_cell", "monocyte", "nk_cell"]


# ---------------------------------------------------------------------------
# Cached loaders so Streamlit doesn't re-query SQLite on every interaction
# ---------------------------------------------------------------------------
@st.cache_data
def load_summary() -> pd.DataFrame:
    with db.connect() as conn:
        return db.summary_table(conn)


@st.cache_data
def load_responder_frame() -> pd.DataFrame:
    with db.connect() as conn:
        return db.responder_comparison_frame(conn)


@st.cache_data
def load_baseline_samples() -> pd.DataFrame:
    with db.connect() as conn:
        return db.baseline_miraclib_melanoma_pbmc_samples(conn)


@st.cache_data
def load_b_cell_counts(sample_ids: tuple[str, ...]) -> pd.DataFrame:
    with db.connect() as conn:
        return db.b_cell_count_for_samples(conn, list(sample_ids))


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Teiko cell-count explorer", layout="wide")
st.title("Teiko cell-count explorer")
st.caption(
    "Interactive view of the analyses from `src/analysis.py`. "
    "All data comes from `teiko.db` (run `python load_data.py` first if missing)."
)

if not (ROOT / "teiko.db").exists():
    st.error("teiko.db not found. Run `python load_data.py` from the repo root first.")
    st.stop()

tab_overview, tab_compare, tab_subset = st.tabs(
    ["Part 2 — Per-sample frequencies", "Part 3 — Responders vs non-responders", "Part 4 — Baseline subset"]
)

# ---------------------------------------------------------------------------
# Tab 1: Part 2 summary table
# ---------------------------------------------------------------------------
with tab_overview:
    st.subheader("Relative frequency of each cell population per sample")
    summary = load_summary()

    col_filt, col_show = st.columns([1, 3])
    with col_filt:
        pops = st.multiselect("Populations", POPULATIONS, default=POPULATIONS)
        sample_search = st.text_input("Filter sample id contains")
    filtered = summary[summary["population"].isin(pops)]
    if sample_search:
        filtered = filtered[filtered["sample"].str.contains(sample_search, case=False, na=False)]
    with col_show:
        st.metric("Rows shown", f"{len(filtered):,}")
        st.dataframe(filtered, use_container_width=True, hide_index=True, height=500)
        st.download_button(
            "Download as CSV",
            filtered.to_csv(index=False).encode(),
            file_name="summary_table.csv",
            mime="text/csv",
        )


# ---------------------------------------------------------------------------
# Tab 2: Part 3 responder comparison
# ---------------------------------------------------------------------------
with tab_compare:
    st.subheader("Melanoma + miraclib + PBMC: responders (yes) vs non-responders (no)")
    long_df = load_responder_frame()

    # Stats panel
    from scipy import stats

    rows = []
    for pop in POPULATIONS:
        sub = long_df[long_df["population"] == pop]
        r = sub.loc[sub["response"] == "yes", "percentage"].to_numpy()
        nr = sub.loc[sub["response"] == "no", "percentage"].to_numpy()
        u, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
        rows.append({"population": pop, "n_resp": len(r), "n_non_resp": len(nr),
                     "median_resp_%": round(float(pd.Series(r).median()), 3),
                     "median_non_resp_%": round(float(pd.Series(nr).median()), 3),
                     "p_value": p})
    stats_df = pd.DataFrame(rows)
    stats_df["p_value_bonferroni"] = (stats_df["p_value"] * len(stats_df)).clip(upper=1.0)
    stats_df["significant (alpha=0.05)"] = stats_df["p_value_bonferroni"] < 0.05

    col_plot, col_stats = st.columns([3, 2])
    with col_plot:
        fig = px.box(
            long_df,
            x="population",
            y="percentage",
            color="response",
            category_orders={"population": POPULATIONS, "response": ["yes", "no"]},
            points="all",
            title="Cell population relative frequency by response",
        )
        fig.update_layout(yaxis_title="Relative frequency (%)", xaxis_title="Population")
        st.plotly_chart(fig, use_container_width=True)
    with col_stats:
        st.markdown("**Mann-Whitney U test, Bonferroni-corrected across 5 populations**")
        st.dataframe(
            stats_df.style.format({"p_value": "{:.4f}", "p_value_bonferroni": "{:.4f}"}),
            use_container_width=True,
            hide_index=True,
        )
        st.caption(
            "Mann-Whitney U is used because cell-frequency distributions are typically "
            "skewed and small-sample, violating t-test normality assumptions. "
            "Bonferroni correction guards against false positives across the 5 simultaneous tests."
        )

# ---------------------------------------------------------------------------
# Tab 3: Part 4 baseline subset
# ---------------------------------------------------------------------------
with tab_subset:
    st.subheader("Baseline cohort: melanoma + miraclib + PBMC, time_from_treatment_start = 0")
    samples = load_baseline_samples()
    subject_attrs = samples.drop_duplicates(subset="subject_id")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Total samples", len(samples))
    c2.metric("Distinct subjects", subject_attrs["subject_id"].nunique())
    c3.metric("Responders", int((subject_attrs["response"] == "yes").sum()))
    c4.metric("Non-responders", int((subject_attrs["response"] == "no").sum()))

    col_a, col_b, col_c = st.columns(3)
    with col_a:
        st.markdown("**Samples per project**")
        st.dataframe(samples["project_id"].value_counts().rename_axis("project").reset_index(name="samples"),
                     hide_index=True, use_container_width=True)
    with col_b:
        st.markdown("**Subjects by response**")
        st.dataframe(subject_attrs["response"].value_counts(dropna=False).rename_axis("response").reset_index(name="subjects"),
                     hide_index=True, use_container_width=True)
    with col_c:
        st.markdown("**Subjects by sex**")
        st.dataframe(subject_attrs["sex"].value_counts().rename_axis("sex").reset_index(name="subjects"),
                     hide_index=True, use_container_width=True)

    st.markdown("**Average B-cell count for melanoma male responders at baseline**")
    male_resp = samples[(samples["sex"] == "M") & (samples["response"] == "yes")]
    bcells = load_b_cell_counts(tuple(male_resp["sample_id"].tolist()))
    avg_bcells = float(bcells["count"].mean()) if len(bcells) else float("nan")
    st.metric("Average B-cell count", f"{avg_bcells:.2f}", help=f"Computed over {len(male_resp)} samples")

    with st.expander("Show baseline samples table"):
        st.dataframe(samples, hide_index=True, use_container_width=True)
