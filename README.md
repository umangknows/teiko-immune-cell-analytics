# Teiko Immune Cell Analytics

This project loads clinical-trial immune cell counts into SQLite, computes per-sample cell population frequencies, compares melanoma PBMC miraclib responders against non-responders, and serves an interactive Streamlit dashboard.

## Run

```bash
make setup
make pipeline
make dashboard
```

The local dashboard runs at `http://localhost:8501`. The app is ready for free deployment on Streamlit Community Cloud from a GitHub repository using `dashboard.py` as the entry point.

## Outputs

`make pipeline` creates:

- `teiko.db`
- `outputs/frequency_summary.csv`
- `outputs/responder_stats.csv`
- `outputs/candidate_signal_scores.csv`
- `outputs/exploratory_prediction_summary.csv`
- `outputs/exploratory_prediction_importance.csv`
- `outputs/analysis_notes.md`
- `outputs/responder_boxplot_baseline.png`
- `outputs/responder_boxplot_all_timepoints.png`
- `outputs/baseline_subset_summary.xlsx`
- `outputs/melanoma_male_bcell_answer.txt`

## Schema

The SQLite database uses a compact normalized schema:

- `projects`: one row per project.
- `subjects`: subject-level metadata including condition, treatment, response, age, and sex.
- `samples`: sample-level metadata including sample type and treatment timepoint.
- `cell_populations`: one row per immune cell population.
- `cell_counts`: long-form cell counts by sample and population.

This avoids hardcoding the five provided cell types into the database model. If future studies add more projects, subjects, sample types, timepoints, or immune populations, the database grows by adding rows rather than changing table structure. The `sample_cell_frequencies` view exposes analysis-ready rows with sample metadata, total count, population count, and relative percentage.

## Analysis Design

Part 2 uses the `sample_cell_frequencies` view and writes the exact requested columns: `sample`, `total_count`, `population`, `count`, and `percentage`.

Part 3 has two layers. The primary inferential comparison uses baseline melanoma PBMC miraclib samples so each subject contributes one independent sample. An exploratory all-timepoint comparison is also reported because the task references trial progression, but repeated measures should not be treated as equally strong inferential evidence without a longitudinal model.

For each population, the pipeline reports responder/non-responder medians, median differences, Mann-Whitney U p-values, Benjamini-Hochberg adjusted p-values, and rank-biserial effect sizes. It also trains a small exploratory logistic-regression model on baseline relative frequencies to provide a candidate predictive signal, reported with cross-validated AUC and coefficient importance. This is intentionally framed as exploratory, not clinically validated.

Part 4 writes baseline subset counts by project, response, and sex, plus the requested average B-cell count for melanoma male responders at time 0 rounded to two decimals.

## Dashboard Design

The dashboard is organized for a project review rather than raw data inspection:

- `Overview`: trial composition, cohort balance, coverage, and the main statistical readout.
- `Cell Frequencies`: filterable frequency trends, sample-type comparisons, and treatment heatmaps.
- `Miraclib Response`: effect sizes, adjusted p-values, boxplots, and exploratory prediction importance.
- `Subset Queries`: visual answers to the baseline melanoma PBMC miraclib subset questions, with audit tables available on demand.

Raw tables are kept in expandable audit sections so the first view emphasizes insight while still preserving traceability.

## Code Structure

- `load_data.py`: required root-level database initializer and CSV loader.
- `run_pipeline.py`: reproducible end-to-end analysis runner.
- `dashboard.py`: Streamlit dashboard.
- `teiko/validation.py`: input checks for required columns, uniqueness, numeric counts, and invalid totals.
- `teiko/db.py`: database connection and loading helpers.
- `teiko/analysis.py`: SQL-backed analysis queries.
- `teiko/stats.py`: responder statistical comparison.
- `teiko/prediction.py`: exploratory response-signal model.
- `teiko/plots.py`: static plot generation.

## Deployment

Recommended free deployment: Streamlit Community Cloud.

1. Push this repository to GitHub.
2. Go to `https://share.streamlit.io`.
3. Create an app from the repository.
4. Use `dashboard.py` as the entry point.
5. Add the resulting `streamlit.app` URL here.
