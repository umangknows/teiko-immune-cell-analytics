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
POPULATION_COLORS = ["#E45756", "#F2CF5B", "#54A24B", "#4C78A8", "#72B7B2"]


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


def scoped_signal_summary(summary: pd.DataFrame, analysis: str) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame()
    out = summary[summary["analysis"] == analysis].copy()
    out["significance"] = out["significant_fdr_0_05"].map(
        {True: "FDR significant", False: "Not significant"}
    )
    return out.sort_values("effect_size_rank_biserial", key=lambda s: s.abs(), ascending=False)


def short_direction(direction: str) -> str:
    if direction == "Higher in responders":
        return "R higher"
    return "NR higher"


def population_order_by_median(df: pd.DataFrame) -> list[str]:
    return (
        df.groupby("population")["percentage"]
        .median()
        .sort_values(ascending=False)
        .index.tolist()
    )


def population_color_map(population_order: list[str]) -> dict[str, str]:
    return {
        population: POPULATION_COLORS[index % len(POPULATION_COLORS)]
        for index, population in enumerate(population_order)
    }


def frequency_takeaway(filtered: pd.DataFrame, median_by_time: pd.DataFrame) -> str:
    population_order = population_order_by_median(filtered)
    dominant = population_order[0]
    changes = []
    for population, group in median_by_time.groupby("population"):
        ordered = group.sort_values("time_from_treatment_start")
        first = float(ordered.iloc[0]["percentage"])
        last = float(ordered.iloc[-1]["percentage"])
        changes.append((population, last - first, first, last))
    largest = max(changes, key=lambda item: abs(item[1]))
    direction = "increased" if largest[1] > 0 else "decreased"
    return (
        f"In simple terms, median relative frequency is the typical share of a sample made up by a given cell population. "
        f"In the selected cohort, {dominant} has the highest median relative frequency, meaning it is usually the largest part of the measured immune-cell mix. "
        f"The largest day 0 to day 14 movement is {largest[0]}, which {direction} by "
        f"{abs(largest[1]):.2f} percentage points ({largest[2]:.2f}% to {largest[3]:.2f}%). "
        "Small movements suggest composition is fairly stable under the current filters."
    )


def response_scope_note(baseline_only: bool) -> str:
    if baseline_only:
        return (
            "Scope: melanoma patients receiving miraclib, PBMC samples only, day 0 baseline. Each subject contributes one sample, "
            "so this is the cleanest responder/non-responder comparison."
        )
    return (
        "Scope: melanoma patients receiving miraclib, PBMC samples only, pooling days 0, 7, and 14. That makes the view useful for "
        "pattern discovery, but the same subject can contribute multiple rows."
    )


def median_change_from_baseline(responder: pd.DataFrame) -> pd.DataFrame:
    medians = (
        responder.groupby(
            ["time_from_treatment_start", "population", "response_label"],
            as_index=False,
        )["percentage"]
        .median()
        .sort_values(["population", "response_label", "time_from_treatment_start"])
    )
    baseline = medians[medians["time_from_treatment_start"] == 0][
        ["population", "response_label", "percentage"]
    ].rename(columns={"percentage": "baseline_percentage"})
    out = medians.merge(baseline, on=["population", "response_label"], how="left")
    out["change_from_baseline"] = out["percentage"] - out["baseline_percentage"]
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
    compare_by = st.selectbox(
        "Compare trend by",
        ["Overall", "Treatment", "Condition", "Sample type", "Response"],
        index=1,
        help="Use this when multiple treatments, conditions, sample types, or response groups are selected. Overall collapses selected groups into one median trend.",
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
        population_order = population_order_by_median(filtered)
        color_map = population_color_map(population_order)
        median_by_time = (
            filtered.groupby(["population", "time_from_treatment_start"], as_index=False)[
                "percentage"
            ]
            .median()
            .sort_values(["population", "time_from_treatment_start"])
        )
        compare_columns = {
            "Treatment": "treatment",
            "Condition": "condition",
            "Sample type": "sample_type",
            "Response": "response_label",
        }
        if compare_by == "Overall":
            trend_fig = px.line(
                median_by_time,
                x="time_from_treatment_start",
                y="percentage",
                color="population",
                markers=True,
                title="How does median immune composition change over treatment time?",
                category_orders={"population": population_order},
                color_discrete_map=color_map,
                labels={
                    "percentage": "Median relative frequency (%)",
                    "time_from_treatment_start": "Days from treatment start",
                },
            )
        else:
            compare_col = compare_columns[compare_by]
            trend_by_group = (
                filtered.groupby(
                    ["population", "time_from_treatment_start", compare_col],
                    as_index=False,
                )["percentage"]
                .median()
                .sort_values([compare_col, "population", "time_from_treatment_start"])
            )
            trend_fig = px.line(
                trend_by_group,
                x="time_from_treatment_start",
                y="percentage",
                color=compare_col,
                facet_col="population",
                facet_col_wrap=3,
                markers=True,
                title=f"How does median immune composition change over time by {compare_by.lower()}?",
                category_orders={"population": population_order},
                labels={
                    "percentage": "Median relative frequency (%)",
                    "time_from_treatment_start": "Days from treatment start",
                    compare_col: compare_by,
                },
            )
            trend_fig.update_yaxes(matches=None)
        st.plotly_chart(style_figure(trend_fig), width="stretch")
        st.info(frequency_takeaway(filtered, median_by_time))
        st.caption(
            "If multiple groups are selected, use `Compare trend by` to avoid mixing them into one line. "
            "For example, compare by Treatment to see miraclib and phauximab separately."
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
                    category_orders={"population": population_order},
                    color_discrete_map=color_map,
                )
            ),
            width="stretch",
        )
        left.caption(
            "This ranks the selected samples by median population percentage. It is a simple composition summary, not a significance test."
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
        right.caption(
            "This compares median composition between PBMC and whole-blood samples when both are present in the filtered data."
        )

    with st.expander("Audit table: exact Part 2 output columns"):
        st.dataframe(load_frequency_summary().head(500), width="stretch", hide_index=True)

with response_tab:
    st.subheader("Do miraclib responders differ from non-responders?")
    st.caption(
        "Scope for every chart on this tab: melanoma patients receiving miraclib, PBMC samples only, responders vs non-responders."
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
    signal_summary = scoped_signal_summary(
        load_csv_output("response_signal_summary.csv"), analysis_key
    )
    change_summary = with_response_labels(load_csv_output("change_from_baseline_summary.csv"))
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
        st.caption(response_scope_note(baseline_only))
        best = stats.sort_values("signal_rank").iloc[0]
        significant_count = int(stats["significant_fdr_0_05"].sum())
        c1, c2, c3, c4 = st.columns(4)
        c1.metric(
            "Top candidate",
            best["population"],
            help="Population with the largest absolute responder/non-responder effect size in this scope.",
        )
        c2.metric(
            "Direction",
            short_direction(best["direction"]),
            help="R higher means responder median is higher. NR higher means non-responder median is higher.",
        )
        c3.metric(
            "Median difference",
            f"{best['median_difference_pct']:.2f} pts",
            help="Responder median percentage minus non-responder median percentage. Negative means non-responders are higher.",
        )
        c4.metric(
            "FDR significant",
            significant_count,
            help="Number of populations with Benjamini-Hochberg adjusted p-value below 0.05.",
        )

        if not signal_summary.empty:
            st.subheader("Response signal summary")
            st.caption(
                "Compact readout for Bob/Yah: median responder frequency, median non-responder frequency, difference, effect size, adjusted p-value, and interpretation."
            )
            signal_fig = px.bar(
                signal_summary.sort_values("median_difference_pct"),
                x="median_difference_pct",
                y="population",
                orientation="h",
                color="interpretation",
                title=f"Miraclib melanoma PBMC: median responder minus non-responder difference ({analysis_scope})",
                labels={
                    "median_difference_pct": "Median difference (percentage points)",
                    "population": "Cell population",
                },
                color_discrete_map={
                    "Higher in responders": "#2B8CBE",
                    "Higher in non-responders": "#F03B20",
                },
                hover_data={
                    "responder_median_pct": ":.3f",
                    "non_responder_median_pct": ":.3f",
                    "adjusted_p_value": ":.4f",
                    "significance": True,
                },
            )
            signal_fig.add_vline(x=0, line_width=1, line_dash="dash", line_color="gray")
            st.plotly_chart(style_figure(signal_fig), width="stretch")
            st.dataframe(
                signal_summary[
                    [
                        "population",
                        "responder_median_pct",
                        "non_responder_median_pct",
                        "median_difference_pct",
                        "effect_size_rank_biserial",
                        "adjusted_p_value",
                        "significance",
                        "interpretation",
                    ]
                ],
                width="stretch",
                hide_index=True,
            )
            st.info(
                "This table is the most compact answer to Part 3. It shows the direction and strength of each candidate signal, while preserving the statistical caution that none are FDR significant in the current dataset."
            )

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
        effect_fig.update_layout(
            title=f"Miraclib melanoma PBMC: effect size direction ({analysis_scope})"
        )
        st.plotly_chart(style_figure(effect_fig), width="stretch")
        st.info(
            "How to read this: this chart is only for melanoma patients treated with miraclib using PBMC samples. "
            "Bars to the right of zero are higher in responders; bars to the left are higher in non-responders. "
            "Longer bars mean larger separation between groups. In this dataset, the bars are small, so the observed differences are weak."
        )

        significance = stats.assign(
            minus_log10_adjusted_p=lambda x: -np.log10(x["adjusted_p_value"])
        )
        threshold = -np.log10(0.05)
        sig_fig = px.bar(
            significance.sort_values("minus_log10_adjusted_p", ascending=False),
            x="population",
            y="minus_log10_adjusted_p",
            color="evidence",
            title=f"Miraclib melanoma PBMC: statistical evidence ({analysis_scope})",
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
        st.info(
            "How to read this: the y-axis is transformed, not a raw p-value. It shows -log10(FDR-adjusted p-value), "
            "so stronger evidence appears taller. The dashed line corresponds to adjusted p=0.05 after correcting for five tested populations. "
            "A bar must rise above that line to be called statistically significant. Current bars stay below the threshold."
        )

    box = px.box(
        responder,
        x="population",
        y="percentage",
        color="response_label",
        points=False,
        title=f"Miraclib melanoma PBMC: responder vs non-responder frequency distributions ({analysis_scope})",
        labels={"percentage": "Relative frequency (%)", "population": "Cell population"},
        color_discrete_map={"Responder": "#2B8CBE", "Non-responder": "#F03B20"},
    )
    st.plotly_chart(style_figure(box), width="stretch")
    st.caption(
        "Each box summarizes the distribution of relative frequencies for one population. Clear vertical separation between responder and non-responder boxes would suggest a stronger response-associated signal."
    )

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
            title="Miraclib melanoma PBMC: median frequency over time by response group",
            labels={
                "percentage": "Median relative frequency (%)",
                "time_from_treatment_start": "Days from treatment start",
            },
            color_discrete_map={"Responder": "#2B8CBE", "Non-responder": "#F03B20"},
        )
        trend_fig.update_yaxes(matches=None)
        st.plotly_chart(style_figure(trend_fig), width="stretch")
        st.info(
            "This view shows median relative frequency over days 0, 7, and 14 separately for responders and non-responders. "
            "It helps reveal whether a signal appears after treatment starts, even if the baseline-only comparison is weak."
        )

        delta = change_summary if not change_summary.empty else median_change_from_baseline(responder)
        if "change_from_baseline_pct_points" in delta.columns:
            delta = delta.rename(columns={"change_from_baseline_pct_points": "change_from_baseline"})
        delta_fig = px.line(
            delta,
            x="time_from_treatment_start",
            y="change_from_baseline",
            color="response_label",
            facet_col="population",
            facet_col_wrap=3,
            markers=True,
            title="Miraclib melanoma PBMC: median change from day 0 baseline",
            labels={
                "change_from_baseline": "Change from baseline (percentage points)",
                "time_from_treatment_start": "Days from treatment start",
            },
            color_discrete_map={"Responder": "#2B8CBE", "Non-responder": "#F03B20"},
        )
        delta_fig.update_yaxes(matches=None)
        delta_fig.add_hline(y=0, line_width=1, line_dash="dash", line_color="gray")
        st.plotly_chart(style_figure(delta_fig), width="stretch")
        st.info(
            "Why this matters: absolute frequencies can look flat even when treatment changes a population relative to its own baseline. "
            "This chart resets each response group to zero at day 0, then shows whether responders and non-responders drift differently after treatment starts."
        )

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
        st.caption(
            "These coefficients show which baseline populations the exploratory logistic model leaned on most. The AUC near 0.50 says that leaning did not translate into useful predictive performance."
        )

    with st.expander("Audit table: statistical results"):
        st.dataframe(stats, width="stretch", hide_index=True)

with query_tab:
    st.subheader("Required Part 4 database query results")
    st.caption(
        "These visuals answer the exact subset questions for baseline melanoma PBMC samples treated with miraclib."
    )
    st.info(
        "This tab is intentionally narrow: Part 4 specifically asks for baseline melanoma PBMC miraclib samples, "
        "counts by project/response/sex, and the average B-cell count for melanoma male responders at time 0."
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
        help="This exact two-decimal value is requested in Part 4. It uses melanoma male responders at time 0 across all sample and treatment types, as specified by the prompt.",
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
