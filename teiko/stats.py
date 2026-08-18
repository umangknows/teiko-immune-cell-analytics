from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def _rank_biserial(u_statistic: float, n_yes: int, n_no: int) -> float:
    return 2.0 * u_statistic / (n_yes * n_no) - 1.0


def compare_responders(freq: pd.DataFrame, analysis_label: str) -> pd.DataFrame:
    rows = []
    for population, group in freq.groupby("population"):
        yes = group.loc[group["response"] == "yes", "percentage"].dropna()
        no = group.loc[group["response"] == "no", "percentage"].dropna()

        if len(yes) < 3 or len(no) < 3:
            rows.append(
                {
                    "analysis": analysis_label,
                    "population": population,
                    "n_responders": len(yes),
                    "n_non_responders": len(no),
                    "median_responder_pct": np.nan,
                    "median_non_responder_pct": np.nan,
                    "median_difference_pct": np.nan,
                    "p_value": np.nan,
                    "effect_size_rank_biserial": np.nan,
                }
            )
            continue

        test = mannwhitneyu(yes, no, alternative="two-sided")
        rows.append(
            {
                "analysis": analysis_label,
                "population": population,
                "n_responders": len(yes),
                "n_non_responders": len(no),
                "median_responder_pct": round(float(yes.median()), 4),
                "median_non_responder_pct": round(float(no.median()), 4),
                "median_difference_pct": round(float(yes.median() - no.median()), 4),
                "p_value": float(test.pvalue),
                "effect_size_rank_biserial": round(_rank_biserial(float(test.statistic), len(yes), len(no)), 4),
            }
        )

    out = pd.DataFrame(rows)
    valid = out["p_value"].notna()
    out["adjusted_p_value"] = np.nan
    out.loc[valid, "adjusted_p_value"] = multipletests(
        out.loc[valid, "p_value"], method="fdr_bh"
    )[1]
    out["adjusted_p_value"] = out["adjusted_p_value"].round(6)
    out["p_value"] = out["p_value"].round(6)
    out["significant_fdr_0_05"] = out["adjusted_p_value"] < 0.05
    out["signal_rank"] = (
        out["effect_size_rank_biserial"].abs().fillna(0).rank(ascending=False, method="first").astype(int)
    )
    return out.sort_values(["analysis", "signal_rank"])


def signal_scores(stats: pd.DataFrame) -> pd.DataFrame:
    scored = stats.copy()
    scored["signal_score"] = (
        scored["effect_size_rank_biserial"].abs().fillna(0) * 100
        + (-np.log10(scored["adjusted_p_value"].clip(lower=1e-300))).fillna(0)
    ).round(3)
    return scored.sort_values("signal_score", ascending=False)

