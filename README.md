# Teiko Immune Cell Analytics

[![Teiko pipeline](https://github.com/umangknows/teiko-immune-cell-analytics/actions/workflows/ci.yml/badge.svg)](https://github.com/umangknows/teiko-immune-cell-analytics/actions/workflows/ci.yml)

This project loads clinical-trial immune cell counts into SQLite, computes per-sample cell population frequencies, compares melanoma PBMC miraclib responders against non-responders, and serves an interactive Streamlit dashboard.

Repository: `https://github.com/umangknows/teiko-immune-cell-analytics`

Dashboard: `https://teiko-immune-cell-analytics.streamlit.app/`

## Run

```bash
make setup
make pipeline
make dashboard
```

The local dashboard runs at `http://localhost:8501`. The public Streamlit dashboard is available at `https://teiko-immune-cell-analytics.streamlit.app/`.

## Reviewer Quickstart

In GitHub Codespaces or a fresh Linux environment:

```bash
make setup
make pipeline
make test
make dashboard
```

Expected sanity checks after `make pipeline`:

- `teiko.db` exists in the repository root.
- `outputs/frequency_summary.csv` has 52,500 rows and the exact columns `sample`, `total_count`, `population`, `count`, `percentage`.
- `outputs/melanoma_male_bcell_answer.txt` contains `10206.15`.
- `outputs/responder_stats.csv` reports responder/non-responder statistics for all five populations.

## Assignment Checklist

| Requirement | Where implemented |
| --- | --- |
| Root-level `load_data.py` initializes SQLite and loads all CSV rows | `load_data.py` |
| Running `python load_data.py` creates a root-level `.db` file | `teiko.db` |
| Relational SQLite schema | `teiko/schema.py` |
| Part 2 frequency table with exact required columns | `outputs/frequency_summary.csv` |
| Part 3 melanoma + miraclib + PBMC responder analysis | `outputs/responder_stats.csv`, dashboard Part 3 |
| Responder vs non-responder boxplots | `outputs/responder_boxplot_baseline.png`, `outputs/responder_boxplot_all_timepoints.png`, dashboard Part 3 |
| Significant populations reported with statistics | `outputs/responder_stats.csv`, `outputs/analysis_notes.md` |
| Response signal summary and longitudinal follow-up | `outputs/response_signal_summary.csv`, `outputs/change_from_baseline_summary.csv` |
| Part 4 baseline subset summaries | `outputs/samples_by_project.csv`, `outputs/subjects_by_response.csv`, `outputs/subjects_by_sex.csv` |
| Average B-cell answer with two decimals | `outputs/melanoma_male_bcell_answer.txt` |
| `make setup` installs dependencies | `Makefile` |
| `make pipeline` runs full data pipeline | `Makefile` |
| `make dashboard` starts dashboard server | `Makefile` |
| Dashboard link included | `https://teiko-immune-cell-analytics.streamlit.app/` |
| CI verifies setup, pipeline, tests, and dashboard syntax | `.github/workflows/ci.yml` |

## Outputs

`make pipeline` creates:

- `teiko.db`
- `outputs/frequency_summary.csv`
- `outputs/responder_stats.csv`
- `outputs/response_signal_summary.csv`
- `outputs/change_from_baseline_summary.csv`
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

## Interpretation

The current five-population frequency panel does not show a strong baseline separation between miraclib responders and non-responders. This should be interpreted as an absence of evidence in the available coarse immune-count features, not proof that no biological response signal exists.

Recommended next investigations:

- Add richer biomarkers such as cytokines, gene expression, additional flow-cytometry markers, tumor mutational burden, or clinical covariates.
- Analyze change from baseline across timepoints, not only absolute frequencies.
- Check whether signal appears within subgroups such as project, age, sex, sample type, or disease context.
- Use subject-level longitudinal models if repeated timepoints become central to inference.

## Dashboard Design

The dashboard is organized for a project review rather than raw data inspection:

- `Overview`: trial composition, sample-type coverage, baseline cohort balance, the main statistical readout, recommended next investigations, and an audit preview.
- `Part 2: Cell Frequencies`: filterable relative-frequency trends, an explicit comparison control for treatment/condition/sample type/response, live cohort-specific interpretation, and composition charts, with the exact required output table available for audit.
- `Part 3: Miraclib Response`: baseline and all-timepoint responder comparisons with plain-language labels, explicit melanoma + miraclib + PBMC scope labels, effect-size direction, statistical evidence, distribution plots, change-from-baseline trends, explanatory chart notes, and exploratory prediction context.
- `Part 4: Required Queries`: visual answers to the baseline melanoma PBMC miraclib subset questions, with audit tables available on demand.

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

The dashboard is deployed on Streamlit Community Cloud:

`https://teiko-immune-cell-analytics.streamlit.app/`

Deployment configuration:

- Repository: `umangknows/teiko-immune-cell-analytics`
- Branch: `main`
- Entry point: `dashboard.py`

Streamlit Community Cloud is a free hosted service, so an app may go cold or require a fresh reboot after inactivity. If the public dashboard shows a sleeping/rebuild message or does not load immediately, open the app from the Streamlit workspace or click `Manage app`, then use the Streamlit `Reboot` option. Streamlit documents this reboot flow in its Community Cloud docs: `https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app/reboot-your-app`.
