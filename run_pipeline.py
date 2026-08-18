from __future__ import annotations

import pandas as pd

from teiko.analysis import (
    baseline_subset_summary,
    frequency_summary,
    melanoma_male_responder_bcell_average,
    responder_frequency_data,
)
from teiko.config import DB_FILE, OUTPUT_DIR
from teiko.plots import responder_boxplot
from teiko.prediction import exploratory_prediction
from teiko.stats import compare_responders, signal_scores


def main() -> None:
    if not DB_FILE.exists():
        raise FileNotFoundError("Database not found. Run python load_data.py first.")

    OUTPUT_DIR.mkdir(exist_ok=True)

    frequencies = frequency_summary()
    frequencies.to_csv(OUTPUT_DIR / "frequency_summary.csv", index=False)

    baseline_freq = responder_frequency_data(baseline_only=True)
    all_timepoint_freq = responder_frequency_data(baseline_only=False)

    stats = pd.concat(
        [
            compare_responders(baseline_freq, "primary_baseline"),
            compare_responders(all_timepoint_freq, "exploratory_all_timepoints"),
        ],
        ignore_index=True,
    )
    stats.to_csv(OUTPUT_DIR / "responder_stats.csv", index=False)
    signal_scores(stats).to_csv(OUTPUT_DIR / "candidate_signal_scores.csv", index=False)

    model_summary, model_importance = exploratory_prediction(baseline_freq)
    model_summary.to_csv(OUTPUT_DIR / "exploratory_prediction_summary.csv", index=False)
    model_importance.to_csv(OUTPUT_DIR / "exploratory_prediction_importance.csv", index=False)

    responder_boxplot(
        baseline_freq,
        OUTPUT_DIR / "responder_boxplot_baseline.png",
        "Baseline melanoma PBMC miraclib: responders vs non-responders",
    )
    responder_boxplot(
        all_timepoint_freq,
        OUTPUT_DIR / "responder_boxplot_all_timepoints.png",
        "All timepoints melanoma PBMC miraclib: exploratory comparison",
    )

    subset = baseline_subset_summary()
    with pd.ExcelWriter(OUTPUT_DIR / "baseline_subset_summary.xlsx") as writer:
        for sheet_name, data in subset.items():
            data.to_excel(writer, sheet_name=sheet_name, index=False)
            data.to_csv(OUTPUT_DIR / f"{sheet_name}.csv", index=False)

    avg_b = melanoma_male_responder_bcell_average()
    (OUTPUT_DIR / "melanoma_male_bcell_answer.txt").write_text(f"{avg_b:.2f}\n", encoding="utf-8")
    significant = stats[
        (stats["analysis"] == "primary_baseline") & (stats["significant_fdr_0_05"])
    ]
    top = stats[stats["analysis"] == "primary_baseline"].sort_values("signal_rank").iloc[0]
    notes = [
        "# Analysis Notes",
        "",
        "Primary inference uses baseline melanoma PBMC samples from miraclib-treated subjects to avoid treating repeated timepoints from the same subject as independent evidence.",
        "",
        f"Strongest baseline candidate signal: {top['population']} "
        f"(median difference {top['median_difference_pct']:.4f} percentage points; "
        f"FDR-adjusted p-value {top['adjusted_p_value']:.6f}; "
        f"rank-biserial effect size {top['effect_size_rank_biserial']:.4f}).",
        "",
        "Significant baseline populations at FDR 0.05: "
        + (", ".join(significant["population"].tolist()) if not significant.empty else "none."),
        "",
        f"Exploratory baseline logistic-regression AUC: {model_summary.loc[0, 'cross_validated_auc']:.4f}. This is reported as an exploratory predictive signal only, not a validated clinical model.",
        "",
        "The task text mentions quintazide, but the dataset contains no quintazide treatment or variable. The analysis therefore uses the treatment values present in cell-count.csv.",
    ]
    (OUTPUT_DIR / "analysis_notes.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"Wrote analysis outputs to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
