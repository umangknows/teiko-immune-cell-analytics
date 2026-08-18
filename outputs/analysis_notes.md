# Analysis Notes

Primary inference uses baseline melanoma PBMC samples from miraclib-treated subjects to avoid treating repeated timepoints from the same subject as independent evidence.

Strongest baseline candidate signal: monocyte (median difference -0.6804 percentage points; FDR-adjusted p-value 0.885328; rank-biserial effect size -0.0564).

Significant baseline populations at FDR 0.05: none.

Exploratory baseline logistic-regression AUC: 0.4890. This is reported as an exploratory predictive signal only, not a validated clinical model.

The task text mentions quintazide, but the dataset contains no quintazide treatment or variable. The analysis therefore uses the treatment values present in cell-count.csv.
