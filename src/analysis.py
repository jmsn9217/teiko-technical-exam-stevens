"""Run all analyses (Parts 2-4) and write outputs to ./outputs/.

Run from the repo root:  python -m src.analysis
"""

from __future__ import annotations

from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")  # headless: works in CI / Codespaces without a display
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

from src import db

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "outputs"

POPULATIONS = ["b_cell", "cd4_t_cell", "cd8_t_cell", "monocyte", "nk_cell"]


# ---------------------------------------------------------------------------
# Part 2: per-sample relative frequency
# ---------------------------------------------------------------------------
def part2_summary_table(conn) -> pd.DataFrame:
    df = db.summary_table(conn)
    out = OUTPUT_DIR / "part2_summary_table.csv"
    df.to_csv(out, index=False)
    print(f"Part 2: wrote {out.relative_to(ROOT)}  ({len(df)} rows)")
    return df


# ---------------------------------------------------------------------------
# Part 3: responder vs non-responder comparison + boxplot + stats
# ---------------------------------------------------------------------------
def part3_responder_comparison(conn) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (long_frame_for_plot, stats_summary)."""
    long_df = db.responder_comparison_frame(conn)

    # Save the per-sample frame so the dashboard / reviewer can re-derive everything
    long_out = OUTPUT_DIR / "part3_responder_frequencies.csv"
    long_df.to_csv(long_out, index=False)

    # Stats: Mann-Whitney U per population, Bonferroni-corrected across the 5 tests
    rows = []
    for pop in POPULATIONS:
        sub = long_df[long_df["population"] == pop]
        responders = sub.loc[sub["response"] == "yes", "percentage"].to_numpy()
        non_responders = sub.loc[sub["response"] == "no", "percentage"].to_numpy()
        u, p = stats.mannwhitneyu(responders, non_responders, alternative="two-sided")
        rows.append(
            {
                "population": pop,
                "n_responders": len(responders),
                "n_non_responders": len(non_responders),
                "median_responders_pct": float(pd.Series(responders).median()),
                "median_non_responders_pct": float(pd.Series(non_responders).median()),
                "mannwhitney_u": float(u),
                "p_value": float(p),
            }
        )
    stats_df = pd.DataFrame(rows)
    # Bonferroni: multiply by number of tests, cap at 1.0
    stats_df["p_value_bonferroni"] = (stats_df["p_value"] * len(stats_df)).clip(upper=1.0)
    stats_df["significant_alpha_0.05"] = stats_df["p_value_bonferroni"] < 0.05

    stats_out = OUTPUT_DIR / "part3_stats.csv"
    stats_df.to_csv(stats_out, index=False)
    print(f"Part 3: wrote {long_out.relative_to(ROOT)}  ({len(long_df)} rows)")
    print(f"Part 3: wrote {stats_out.relative_to(ROOT)}")
    print(stats_df.to_string(index=False))

    # Boxplot
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=long_df,
        x="population",
        y="percentage",
        hue="response",
        hue_order=["yes", "no"],
        order=POPULATIONS,
        ax=ax,
    )
    sns.stripplot(
        data=long_df,
        x="population",
        y="percentage",
        hue="response",
        hue_order=["yes", "no"],
        order=POPULATIONS,
        dodge=True,
        size=3,
        alpha=0.5,
        ax=ax,
        legend=False,
    )
    # Annotate significant populations with an asterisk
    for i, pop in enumerate(POPULATIONS):
        if bool(stats_df.loc[stats_df["population"] == pop, "significant_alpha_0.05"].iloc[0]):
            top = long_df.loc[long_df["population"] == pop, "percentage"].max()
            ax.text(i, top + 1.5, "*", ha="center", va="bottom", fontsize=18, fontweight="bold")

    ax.set_title("Cell population relative frequency: melanoma + miraclib + PBMC\nresponders (yes) vs non-responders (no)")
    ax.set_xlabel("Cell population")
    ax.set_ylabel("Relative frequency (%)")
    ax.legend(title="Response")
    fig.tight_layout()
    plot_out = OUTPUT_DIR / "part3_boxplot.png"
    fig.savefig(plot_out, dpi=150)
    plt.close(fig)
    print(f"Part 3: wrote {plot_out.relative_to(ROOT)}")
    return long_df, stats_df


# ---------------------------------------------------------------------------
# Part 4: subset analysis on baseline melanoma + miraclib + PBMC samples
# ---------------------------------------------------------------------------
def part4_baseline_subset(conn) -> dict:
    samples = db.baseline_miraclib_melanoma_pbmc_samples(conn)

    by_project = samples["project_id"].value_counts().sort_index().to_dict()
    # Subject-level counts: collapse duplicates if a subject ever has >1 baseline sample
    subject_attrs = samples.drop_duplicates(subset="subject_id")
    by_response = subject_attrs.dropna(subset=["response"])["response"].value_counts().to_dict()
    by_sex = subject_attrs["sex"].value_counts().to_dict()

    # Average B cells (raw count) for melanoma MALE responders at time=0
    male_resp_samples = samples[(samples["sex"] == "M") & (samples["response"] == "yes")]
    bcells = db.b_cell_count_for_samples(conn, male_resp_samples["sample_id"].tolist())
    avg_bcell_male_responders = float(bcells["count"].mean()) if len(bcells) else float("nan")

    samples_out = OUTPUT_DIR / "part4_baseline_samples.csv"
    samples.to_csv(samples_out, index=False)

    summary = {
        "n_samples_total": int(len(samples)),
        "samples_per_project": {str(k): int(v) for k, v in by_project.items()},
        "subjects_by_response": {str(k): int(v) for k, v in by_response.items()},
        "subjects_by_sex": {str(k): int(v) for k, v in by_sex.items()},
        "avg_b_cell_count_male_responders_baseline": round(avg_bcell_male_responders, 2),
        "n_male_responder_baseline_samples": int(len(male_resp_samples)),
    }

    summary_out = OUTPUT_DIR / "part4_summary.json"
    summary_out.write_text(json.dumps(summary, indent=2))
    print(f"Part 4: wrote {samples_out.relative_to(ROOT)}  ({len(samples)} rows)")
    print(f"Part 4: wrote {summary_out.relative_to(ROOT)}")
    print(json.dumps(summary, indent=2))
    return summary


def main() -> int:
    OUTPUT_DIR.mkdir(exist_ok=True)
    sns.set_theme(style="whitegrid")
    with db.connect() as conn:
        part2_summary_table(conn)
        print()
        part3_responder_comparison(conn)
        print()
        part4_baseline_subset(conn)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
