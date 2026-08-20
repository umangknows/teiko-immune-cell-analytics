from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from teiko.analysis import (
    baseline_subset_summary,
    dataset_overview,
    frequency_summary,
    frequency_summary_with_metadata,
    melanoma_male_responder_bcell_average,
    responder_frequency_data,
    sample_metadata,
)
from teiko.config import DB_FILE, OUTPUT_DIR


st.set_page_config(page_title="Teiko Immune Cell Analytics", layout="wide")

RESPONSE_LABELS = {"yes": "Responder", "no": "Non-responder"}


def with_response_labels(df: pd.DataFrame) -> pd.DataFrame:
    if "response" not in df.columns:
        return df
    out = df.copy()
    out["response_label"] = out["response"].map(RESPONSE_LABELS).fillna("Not applicable")
    return out


@st.cache_data
def load_frequency_summary() -> pd.DataFrame:
    return frequency_summary()


@st.cache_data
def load_full_frequencies() -> pd.DataFrame:
    return with_response_labels(frequency_summary_with_metadata())


@st.cache_data
def load_samples() -> pd.DataFrame:
    return with_response_labels(sample_metadata())


@st.cache_data
def load_responder_data(baseline_only: bool) -> pd.DataFrame:
    return with_response_labels(responder_frequency_data(baseline_only=baseline_only))


@st.cache_data
def load_csv_output(name: str) -> pd.DataFrame:
    path = OUTPUT_DIR / name
    if path.exists():
        return pd.read_csv(path)
    return pd.DataFrame()


def style_figure(fig: go.Figure) -> go.Figure:
    fig.update_layout(
        margin=dict(l=10, r=10, t=54, b=28),
        legend_title_text="",
        font=dict(size=13),
    )
    return fig


def scoped_stats(stats: pd.DataFrame, analysis: str) -> pd.DataFrame:
    if stats.empty:
        return pd.DataFrame()
    out = stats[stats["analysis"] == analysis].copy()
    out["direction"] = out["effect_size_rank_biserial"].apply(
        lambda value: "Higher in responders" if value > 0 else "Higher in non-responders"
    )
    out["evidence"] = out["significant_fdr_0_05"].map(
        {True: "FDR significant", False: "Not significant"}
    )
    return out


st.title("Teiko Immune Cell Analytics")
st.caption(
    "A reproducible clinical-trial data pipeline and dashboard for immune cell population analysis."
)

if not DB_FILE.exists():
    st.error("Database not found. Run `make pipeline` before launching the dashboard.")
    st.stop()

overview_tab, frequency_tab, response_tab, query_tab = st.tabs(
    [
        "Overview",
        "Part 2: Cell Frequencies",
        "Part 3: Miraclib Response",
        "Part 4: Required Queries",
    ]
)

with overview_tab:
    overview = dataset_overview()
    samples = load_samples()
    baseline_freq = load_responder_data(True)
    stats = load_csv_output("responder_stats.csv")
    model_summary = load_csv_output("exploratory_prediction_summary.csv")
    primary = scoped_stats(stats, "primary_baseline")

    cols = st.columns(5)
    cols[0].metric("Projects", f"{overview['projects']:,}")
    cols[1].metric("Subjects", f"{overview['subjects']:,}")
    cols[2].metric("Samples", f"{overview['samples']:,}")
    cols[3].metric("Cell-count rows", f"{overview['cell_count_records']:,}")
    cols[4].metric("Primary cohort", f"{baseline_freq['sample'].nunique():,} samples")

    if not primary.empty:
        top = primary.sort_values("signal_rank").iloc[0]
        auc = model_summary.loc[0, "cross_validated_auc"] if not model_summary.empty else "N/A"
        st.info(
            "Main readout: baseline melanoma PBMC samples from miraclib-treated patients did not show "
            "FDR-significant responder/non-responder differences across the five measured immune populations. "
            f"The strongest baseline candidate was {top['population']} "
            f"({top['direction'].lower()}, median difference {top['median_difference_pct']:.2f} percentage points, "
            f"adjusted p={top['adjusted_p_value']:.3f}). Exploratory baseline prediction AUC was {auc}."
        )

    left, right = st.columns(2)
    composition = samples.groupby(["condition", "treatment"], as_index=False)["sample"].nunique()
    left.plotly_chart(
        style_figure(
            px.bar(
                composition,
                x="condition",
                y="sample",
                color="treatment",
                barmode="stack",
                title="What data do we have?",
                labels={"sample": "Samples", "condition": "Disease/condition"},
                color_discrete_sequence=px.colors.qualitative.Set2,
            )
        ),
        width="stretch",
    )

    sample_types = samples.groupby(["project", "sample_type"], as_index=False)["sample"].nunique()
    right.plotly_chart(
        style_figure(
            px.bar(
                sample_types,
                x="project",
                y="sample",
                color="sample_type",
                barmode="group",
                title="Where are sample types represented?",
                labels={"sample": "Samples", "project": "Project"},
                color_discrete_sequence=["#4C78A8", "#F58518"],
            )
        ),
        width="stretch",
    )

    focus = samples[
        (samples["condition"] == "melanoma")
        & (samples["treatment"] == "miraclib")
        & (samples["sample_type"] == "PBMC")
        & (samples["response"].isin(["yes", "no"]))
    ]
    baseline_focus = focus[focus["time_from_treatment_start"] == 0]
    response_counts = (
        baseline_focus.groupby(["response_label"], as_index=False)["sample"]
        .nunique()
        .rename(columns={"sample": "sample_count"})
    )
    left, right = st.columns([1, 2])
    left.plotly_chart(
        style_figure(
            px.bar(
                response_counts,
                x="response_label",
                y="sample_count",
                color="response_label",
                title="Primary baseline cohort balance",
                labels={"sample_count": "Samples", "response_label": "Response"},
                color_discrete_map={"Responder": "#2B8CBE", "Non-responder": "#F03B20"},
            )
        ),
        width="stretch",
    )
    with right:
        st.subheader("What the current result means")
        st.write(
            "The current five-population frequency panel does not show a strong baseline separation "
            "between miraclib responders and non-responders. That does not prove there is no biology; "
            "it means this dataset may be too coarse to expose it alone."
        )
        st.markdown(
            """
            Recommended next investigations:

            - Add richer biomarkers such as cytokines, gene expression, flow markers, tumor mutational burden, or clinical covariates.
            - Model longitudinal change from baseline instead of only absolute frequencies at a single timepoint.
            - Check whether signals emerge in specific subgroups, projects, ages, sex groups, or sample types.
            - Use subject-level longitudinal models if repeated timepoints become central to the analysis.
            """
        )

    with st.expander("Audit preview: sample metadata"):
        st.dataframe(samples.head(250), width="stretch", hide_index=True)

with frequency_tab:
    st.subheader("Relative frequency of each cell population in each sample")
    st.caption(
        "Part 2 asks for each sample's total count and the percentage contribution of each immune population."
    )

    freq = load_full_frequencies()
    filters = st.columns(5)
    condition = filters[0].multiselect("Condition", sorted(freq["condition"].unique()))
    treatment = filters[1].multiselect("Treatment", sorted(freq["treatment"].unique()))
    sample_type = filters[2].multiselect("Sample type", sorted(freq["sample_type"].unique()))
    timepoint = filters[3].multiselect(
        "Timepoint", sorted(freq["time_from_treatment_start"].unique())
    )
    response = filters[4].multiselect(
        "Response", ["Responder", "Non-responder", "Not applicable"]
    )
    populations = st.multiselect(
        "Population",
        sorted(freq["population"].unique()),
        default=sorted(freq["population"].unique()),
    )

    filtered = freq[freq["population"].isin(populations)]
    for column, values in {
        "condition": condition,
        "treatment": treatment,
        "sample_type": sample_type,
        "time_from_treatment_start": timepoint,
        "response_label": response,
    }.items():
        if values:
            filtered = filtered[filtered[column].isin(values)]

    if filtered.empty:
        st.warning("No rows match the selected filters.")
    else:
        median_by_time = (
            filtered.groupby(["population", "time_from_treatment_start"], as_index=False)[
                "percentage"
            ]
            .median()
            .sort_values(["population", "time_from_treatment_start"])
        )
        st.plotly_chart(
            style_figure(
                px.line(
                    median_by_time,
                    x="time_from_treatment_start",
                    y="percentage",
                    color="population",
                    markers=True,
                    title="How does median immune composition change over treatment time?",
                    labels={
                        "percentage": "Median relative frequency (%)",
                        "time_from_treatment_start": "Days from treatment start",
                    },
                    color_discrete_sequence=px.colors.qualitative.Safe,
                )
            ),
            width="stretch",
        )

        left, right = st.columns(2)
        population_mix = (
            filtered.groupby(["population"], as_index=False)["percentage"]
            .median()
            .sort_values("percentage", ascending=False)
        )
        left.plotly_chart(
            style_figure(
                px.bar(
                    population_mix,
                    x="population",
                    y="percentage",
                    title="Which populations dominate the selected cohort?",
                    labels={"percentage": "Median relative frequency (%)"},
                    color="population",
                    color_discrete_sequence=px.colors.qualitative.Safe,
                )
            ),
            width="stretch",
        )

        by_sample_type = (
            filtered.groupby(["population", "sample_type"], as_index=False)["percentage"].median()
        )
        right.plotly_chart(
            style_figure(
                px.bar(
                    by_sample_type,
                    x="population",
                    y="percentage",
                    color="sample_type",
                    barmode="group",
                    title="Do PBMC and whole-blood samples differ?",
                    labels={"percentage": "Median relative frequency (%)"},
                    color_discrete_sequence=["#4C78A8", "#F58518"],
                )
            ),
            width="stretch",
        )

    with st.expander("Audit table: exact Part 2 output columns"):
        st.dataframe(load_frequency_summary().head(500), width="stretch", hide_index=True)

with response_tab:
    st.subheader("Do miraclib responders differ from non-responders?")
    st.caption(
        "The required comparison is limited to melanoma patients receiving miraclib, using PBMC samples."
    )

    analysis_scope = st.radio(
        "Analysis scope",
        ["Primary: baseline only", "Exploratory: all timepoints"],
        horizontal=True,
    )
    baseline_only = analysis_scope.startswith("Primary")
    analysis_key = "primary_baseline" if baseline_only else "exploratory_all_timepoints"

    responder = load_responder_data(baseline_only)
    stats = scoped_stats(load_csv_output("responder_stats.csv"), analysis_key)
    model_summary = load_csv_output("exploratory_prediction_summary.csv")
    model_importance = load_csv_output("exploratory_prediction_importance.csv")

    if baseline_only:
        st.info(
            "Primary baseline analysis uses time 0 only so each subject contributes one sample. "
            "This is the cleaner view for asking whether baseline cell composition may help predict response."
        )
    else:
        st.warning(
            "Exploratory all-timepoint analysis includes days 0, 7, and 14. It is useful for pattern finding, "
            "but repeated samples from the same subject make it less definitive than the baseline-only test."
        )

    if not stats.empty:
        best = stats.sort_values("signal_rank").iloc[0]
        significant_count = int(stats["significant_fdr_0_05"].sum())
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Top candidate", best["population"])
        c2.metric("Direction", best["direction"])
        c3.metric("Median difference", f"{best['median_difference_pct']:.2f} pts")
        c4.metric("FDR-significant populations", significant_count)

        evidence = stats.sort_values("effect_size_rank_biserial")
        effect_fig = px.bar(
            evidence,
            x="effect_size_rank_biserial",
            y="population",
            orientation="h",
            color="direction",
            title="Effect size direction: which group is higher?",
            labels={
                "effect_size_rank_biserial": "Effect size",
                "population": "Cell population",
            },
            color_discrete_map={
                "Higher in responders": "#2B8CBE",
                "Higher in non-responders": "#F03B20",
            },
            hover_data={
                "median_difference_pct": ":.3f",
                "adjusted_p_value": ":.4f",
                "direction": True,
            },
        )
        effect_fig.add_vline(x=0, line_width=1, line_dash="dash", line_color="gray")
        st.plotly_chart(style_figure(effect_fig), width="stretch")

        significance = stats.assign(
            minus_log10_adjusted_p=lambda x: -np.log10(x["adjusted_p_value"])
        )
        threshold = -np.log10(0.05)
        sig_fig = px.bar(
            significance.sort_values("minus_log10_adjusted_p", ascending=False),
            x="population",
            y="minus_log10_adjusted_p",
            color="evidence",
            title="Statistical evidence after multiple-testing correction",
            labels={
                "minus_log10_adjusted_p": "-log10(FDR-adjusted p-value)",
                "population": "Cell population",
            },
            color_discrete_map={
                "FDR significant": "#2B8CBE",
                "Not significant": "#9AA0A6",
            },
        )
        sig_fig.add_hline(
            y=threshold,
            line_width=1,
            line_dash="dash",
            line_color="#D62728",
            annotation_text="FDR 0.05 threshold",
        )
        st.plotly_chart(style_figure(sig_fig), width="stretch")

    box = px.box(
        responder,
        x="population",
        y="percentage",
        color="response_label",
        points=False,
        title="Responder vs non-responder frequency distributions",
        labels={"percentage": "Relative frequency (%)", "population": "Cell population"},
        color_discrete_map={"Responder": "#2B8CBE", "Non-responder": "#F03B20"},
    )
    st.plotly_chart(style_figure(box), width="stretch")

    if not baseline_only:
        trend = (
            responder.groupby(
                ["time_from_treatment_start", "population", "response_label"], as_index=False
            )["percentage"]
            .median()
            .sort_values(["population", "response_label", "time_from_treatment_start"])
        )
        trend_fig = px.line(
            trend,
            x="time_from_treatment_start",
            y="percentage",
            color="response_label",
            facet_col="population",
            facet_col_wrap=3,
            markers=True,
            title="Median frequency over time by response group",
            labels={
                "percentage": "Median relative frequency (%)",
                "time_from_treatment_start": "Days from treatment start",
            },
            color_discrete_map={"Responder": "#2B8CBE", "Non-responder": "#F03B20"},
        )
        trend_fig.update_yaxes(matches=None)
        st.plotly_chart(style_figure(trend_fig), width="stretch")

    if baseline_only and not model_summary.empty:
        st.subheader("Exploratory predictive signal")
        st.caption(
            "AUC near 0.50 means the five baseline cell-frequency features do not separate responders from non-responders better than chance."
        )
        left, right = st.columns([1, 2])
        left.metric("Cross-validated AUC", f"{model_summary.loc[0, 'cross_validated_auc']:.3f}")
        left.info(model_summary.loc[0, "note"])
        importance_fig = px.bar(
            model_importance.sort_values("absolute_importance"),
            x="absolute_importance",
            y="population",
            orientation="h",
            title="Which baseline features mattered most to the exploratory model?",
            labels={"absolute_importance": "Absolute standardized coefficient"},
            color_discrete_sequence=["#4C78A8"],
        )
        right.plotly_chart(style_figure(importance_fig), width="stretch")

    with st.expander("Audit table: statistical results"):
        st.dataframe(stats, width="stretch", hide_index=True)

with query_tab:
    st.subheader("Required Part 4 database query results")
    st.caption(
        "These visuals answer the exact subset questions for baseline melanoma PBMC samples treated with miraclib."
    )

    subset = baseline_subset_summary()
    avg_b = melanoma_male_responder_bcell_average()
    c1, c2, c3 = st.columns(3)
    c1.plotly_chart(
        style_figure(
            px.bar(
                subset["samples_by_project"],
                x="project",
                y="sample_count",
                title="How many baseline samples from each project?",
                labels={"sample_count": "Samples", "project": "Project"},
                color_discrete_sequence=["#4C78A8"],
            )
        ),
        width="stretch",
    )
    c2.plotly_chart(
        style_figure(
            px.bar(
                subset["subjects_by_response"].assign(
                    response_label=lambda x: x["response"].map(RESPONSE_LABELS)
                ),
                x="response_label",
                y="subject_count",
                title="How many subjects responded?",
                labels={"subject_count": "Subjects", "response_label": "Response"},
                color="response_label",
                color_discrete_map={"Responder": "#2B8CBE", "Non-responder": "#F03B20"},
            )
        ),
        width="stretch",
    )
    c3.plotly_chart(
        style_figure(
            px.bar(
                subset["subjects_by_sex"],
                x="sex",
                y="subject_count",
                title="How many male/female subjects?",
                labels={"subject_count": "Subjects", "sex": "Sex"},
                color_discrete_sequence=["#72B7B2"],
            )
        ),
        width="stretch",
    )

    st.metric(
        "Average B-cell count for melanoma male responders at time 0",
        f"{avg_b:.2f}",
    )

    with st.expander("Audit tables: subset query results"):
        left, middle, right = st.columns(3)
        left.dataframe(subset["samples_by_project"], width="stretch", hide_index=True)
        middle.dataframe(
            subset["subjects_by_response"].assign(
                response_label=lambda x: x["response"].map(RESPONSE_LABELS)
            ),
            width="stretch",
            hide_index=True,
        )
        right.dataframe(subset["subjects_by_sex"], width="stretch", hide_index=True)
