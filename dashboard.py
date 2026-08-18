from __future__ import annotations

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


@st.cache_data
def load_frequency_summary() -> pd.DataFrame:
    return frequency_summary()


@st.cache_data
def load_full_frequencies() -> pd.DataFrame:
    return frequency_summary_with_metadata()


@st.cache_data
def load_samples() -> pd.DataFrame:
    return sample_metadata()


@st.cache_data
def load_responder_data(baseline_only: bool) -> pd.DataFrame:
    return responder_frequency_data(baseline_only=baseline_only)


@st.cache_data
def load_csv_output(name: str) -> pd.DataFrame:
    path = OUTPUT_DIR / name
    if path.exists():
        return pd.read_csv(path)
    return pd.DataFrame()


def style_figure(fig: go.Figure) -> go.Figure:
    fig.update_layout(
        margin=dict(l=10, r=10, t=48, b=20),
        legend_title_text="",
        font=dict(size=13),
    )
    return fig


def responder_delta_heatmap(freq: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        freq.groupby(["time_from_treatment_start", "population", "response"], as_index=False)[
            "percentage"
        ]
        .median()
        .pivot_table(
            index=["time_from_treatment_start", "population"],
            columns="response",
            values="percentage",
        )
        .reset_index()
    )
    grouped["responder_minus_non_responder"] = grouped["yes"] - grouped["no"]
    return grouped.pivot(
        index="population",
        columns="time_from_treatment_start",
        values="responder_minus_non_responder",
    )


st.title("Teiko Immune Cell Analytics")
st.caption("Clinical-trial cell count pipeline, responder analysis, and dashboard.")

if not DB_FILE.exists():
    st.error("Database not found. Run `make pipeline` before launching the dashboard.")
    st.stop()

overview_tab, frequencies_tab, response_tab, subsets_tab = st.tabs(
    ["Overview", "Cell Frequencies", "Miraclib Response", "Subset Queries"]
)

with overview_tab:
    overview = dataset_overview()
    samples = load_samples()
    full_freq = load_full_frequencies()
    baseline_freq = load_responder_data(True)
    stats = load_csv_output("responder_stats.csv")
    model_summary = load_csv_output("exploratory_prediction_summary.csv")
    primary_stats = stats[stats["analysis"] == "primary_baseline"] if not stats.empty else pd.DataFrame()

    cols = st.columns(5)
    for col, (label, value) in zip(cols, overview.items()):
        col.metric(label.replace("_", " ").title(), f"{value:,}")
    cols[-1].metric("Baseline Response Cohort", f"{baseline_freq['sample'].nunique():,}")

    if not primary_stats.empty:
        top = primary_stats.sort_values("signal_rank").iloc[0]
        auc = model_summary.loc[0, "cross_validated_auc"] if not model_summary.empty else "N/A"
        st.info(
            "Primary baseline analysis found no FDR-significant cell-population frequency "
            f"differences. The strongest candidate signal was {top['population']} "
            f"(median difference {top['median_difference_pct']:.2f} pts; adjusted p={top['adjusted_p_value']:.3f}). "
            f"Exploratory baseline prediction AUC: {auc}."
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
                title="Trial composition by condition and treatment",
                labels={"sample": "Samples", "condition": "Condition"},
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
                title="Sample type coverage by project",
                labels={"sample": "Samples", "project": "Project"},
            )
        ),
        width="stretch",
    )

    left, right = st.columns(2)
    response_counts = (
        samples[
            (samples["condition"] == "melanoma")
            & (samples["treatment"] == "miraclib")
            & (samples["sample_type"] == "PBMC")
            & (samples["response"].isin(["yes", "no"]))
        ]
        .groupby(["time_from_treatment_start", "response"], as_index=False)["sample"]
        .nunique()
    )
    left.plotly_chart(
        style_figure(
            px.bar(
                response_counts,
                x="time_from_treatment_start",
                y="sample",
                color="response",
                barmode="group",
                title="Miraclib melanoma PBMC response cohort by timepoint",
                labels={"sample": "Samples", "time_from_treatment_start": "Days from treatment start"},
            )
        ),
        width="stretch",
    )

    median_landscape = (
        full_freq.groupby(["condition", "population"], as_index=False)["percentage"].median()
    )
    right.plotly_chart(
        style_figure(
            px.bar(
                median_landscape,
                x="population",
                y="percentage",
                color="condition",
                barmode="group",
                title="Median immune composition by condition",
                labels={"percentage": "Median relative frequency (%)"},
            )
        ),
        width="stretch",
    )

    with st.expander("Audit preview: frequency summary"):
        st.dataframe(load_frequency_summary().head(100), width="stretch", hide_index=True)

with frequencies_tab:
    freq = load_full_frequencies()
    filters = st.columns(5)
    condition = filters[0].multiselect("Condition", sorted(freq["condition"].unique()))
    treatment = filters[1].multiselect("Treatment", sorted(freq["treatment"].unique()))
    sample_type = filters[2].multiselect("Sample type", sorted(freq["sample_type"].unique()))
    timepoint = filters[3].multiselect(
        "Timepoint", sorted(freq["time_from_treatment_start"].unique())
    )
    response = filters[4].multiselect("Response", sorted(freq["response"].dropna().unique()))
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
        "response": response,
    }.items():
        if values:
            filtered = filtered[filtered[column].isin(values)]

    if not filtered.empty:
        median_by_group = (
            filtered.groupby(["population", "time_from_treatment_start"], as_index=False)[
                "percentage"
            ].median()
        )
        st.plotly_chart(
            style_figure(
                px.line(
                    median_by_group,
                    x="time_from_treatment_start",
                    y="percentage",
                    color="population",
                    markers=True,
                    title="Median relative frequency over treatment time",
                    labels={
                        "percentage": "Median relative frequency (%)",
                        "time_from_treatment_start": "Days from treatment start",
                    },
                )
            ),
            width="stretch",
        )

        left, right = st.columns(2)
        distribution = (
            filtered.groupby(["population", "sample_type"], as_index=False)["percentage"].median()
        )
        left.plotly_chart(
            style_figure(
                px.bar(
                    distribution,
                    x="population",
                    y="percentage",
                    color="sample_type",
                    barmode="group",
                    title="Median frequency by sample type",
                    labels={"percentage": "Median relative frequency (%)"},
                )
            ),
            width="stretch",
        )

        heatmap = filtered.pivot_table(
            index="population",
            columns="treatment",
            values="percentage",
            aggfunc="median",
        )
        right.plotly_chart(
            style_figure(
                px.imshow(
                    heatmap,
                    aspect="auto",
                    color_continuous_scale="Tealrose",
                    title="Median frequency heatmap by treatment",
                    labels={"color": "Median %"},
                )
            ),
            width="stretch",
        )

    with st.expander("Audit table: filtered frequency rows"):
        st.dataframe(filtered, width="stretch", hide_index=True)

with response_tab:
    analysis_scope = st.radio(
        "Analysis scope",
        ["Primary: baseline only", "Exploratory: all timepoints"],
        horizontal=True,
    )
    baseline_only = analysis_scope.startswith("Primary")
    responder = load_responder_data(baseline_only)
    stats = load_csv_output("responder_stats.csv")
    model_summary = load_csv_output("exploratory_prediction_summary.csv")
    model_importance = load_csv_output("exploratory_prediction_importance.csv")
    wanted = "primary_baseline" if baseline_only else "exploratory_all_timepoints"
    scoped_stats = stats[stats["analysis"] == wanted] if not stats.empty else pd.DataFrame()

    if not scoped_stats.empty:
        best = scoped_stats.sort_values("signal_rank").iloc[0]
        c1, c2, c3 = st.columns(3)
        c1.metric("Top candidate population", best["population"])
        c2.metric("Median difference", f"{best['median_difference_pct']:.2f} pts")
        c3.metric("Adjusted p-value", f"{best['adjusted_p_value']:.3f}")

        left, right = st.columns(2)
        effect_fig = px.bar(
            scoped_stats.sort_values("effect_size_rank_biserial"),
            x="effect_size_rank_biserial",
            y="population",
            orientation="h",
            color="effect_size_rank_biserial",
            color_continuous_scale="RdBu",
            title="Effect size by cell population",
            labels={"effect_size_rank_biserial": "Rank-biserial effect size"},
        )
        effect_fig.add_vline(x=0, line_width=1, line_dash="dash", line_color="gray")
        left.plotly_chart(style_figure(effect_fig), width="stretch")

        p_fig = px.scatter(
            scoped_stats,
            x="population",
            y="adjusted_p_value",
            size=scoped_stats["effect_size_rank_biserial"].abs(),
            color="significant_fdr_0_05",
            title="Adjusted p-values with effect-size emphasis",
            labels={"adjusted_p_value": "FDR-adjusted p-value"},
        )
        p_fig.add_hline(y=0.05, line_width=1, line_dash="dash", line_color="gray")
        right.plotly_chart(style_figure(p_fig), width="stretch")

    fig = px.box(
        responder,
        x="population",
        y="percentage",
        color="response",
        points="outliers",
        title="Responder vs non-responder relative frequencies",
    )
    st.plotly_chart(fig, width="stretch")

    if not baseline_only:
        delta = responder_delta_heatmap(responder)
        st.plotly_chart(
            style_figure(
                px.imshow(
                    delta,
                    aspect="auto",
                    color_continuous_scale="RdBu",
                    color_continuous_midpoint=0,
                    title="Median responder minus non-responder frequency by timepoint",
                    labels={"color": "Percentage-point delta"},
                )
            ),
            width="stretch",
        )

    if not model_summary.empty and baseline_only:
        st.subheader("Exploratory Predictive Signal")
        left, right = st.columns([1, 2])
        left.metric("Cross-validated AUC", model_summary.loc[0, "cross_validated_auc"])
        left.caption(model_summary.loc[0, "note"])
        importance_fig = px.bar(
            model_importance.sort_values("absolute_importance"),
            x="absolute_importance",
            y="population",
            orientation="h",
            title="Logistic regression coefficient importance",
            labels={"absolute_importance": "Absolute standardized coefficient"},
        )
        right.plotly_chart(style_figure(importance_fig), width="stretch")

    with st.expander("Audit table: statistical results"):
        st.dataframe(scoped_stats, width="stretch", hide_index=True)

with subsets_tab:
    subset = baseline_subset_summary()
    avg_b = melanoma_male_responder_bcell_average()
    c1, c2, c3 = st.columns(3)
    c1.plotly_chart(
        style_figure(
            px.bar(
                subset["samples_by_project"],
                x="project",
                y="sample_count",
                title="Baseline samples by project",
                labels={"sample_count": "Samples"},
            )
        ),
        width="stretch",
    )
    c2.plotly_chart(
        style_figure(
            px.pie(
                subset["subjects_by_response"],
                names="response",
                values="subject_count",
                title="Subjects by response",
            )
        ),
        width="stretch",
    )
    c3.plotly_chart(
        style_figure(
            px.pie(
                subset["subjects_by_sex"],
                names="sex",
                values="subject_count",
                title="Subjects by sex",
            )
        ),
        width="stretch",
    )
    st.metric(
        "Average B cells for melanoma male responders at time 0",
        f"{avg_b:.2f}",
    )
    with st.expander("Audit tables: subset query results"):
        left, middle, right = st.columns(3)
        left.dataframe(subset["samples_by_project"], width="stretch", hide_index=True)
        middle.dataframe(subset["subjects_by_response"], width="stretch", hide_index=True)
        right.dataframe(subset["subjects_by_sex"], width="stretch", hide_index=True)

