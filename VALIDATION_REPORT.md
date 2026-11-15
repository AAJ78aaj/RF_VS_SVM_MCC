# R-code validation report

The portable script was run in the supplied environment with all four reported spatial geometries enabled. It was checked against the manuscript’s saved numerical analysis records. No manuscript or source input was changed.

| Output | Rows | Maximum absolute numeric difference | Result |
|---|---:|---:|---|
| `archived_seed123_split_counts.csv` | 2 | 0 | PASS |
| `spatial_oof_predictions.csv` | 2042 | 0 | PASS |
| `nested_tuning_results.csv` | 150 | 0 | PASS |
| `spatial_fold_accounting.csv` | 10 | 0 | PASS |
| `spatial_performance_metrics.csv` | 4 | 0 | PASS |
| `spatial_block_bootstrap.csv` | 26 | 0 | PASS |
| `mpr_aspect_rates_with_exact_ci.csv` | 8 | 0 | PASS |
| `mpr_logistic_regression.csv` | 5 | 0 | PASS |
| `mpr_predictor_scaling.csv` | 4 | 0 | PASS |
| `mpr_book_coordinate_overlap.csv` | 1 | 0 | PASS |
| `spatial_geometry_sensitivity.csv` | 8 | 0 | PASS |
| `corrected_original_holdout_metrics.csv` | 4 | 4.44e-16 | PASS |
| `reference_label_moran.csv` | 2 | 1.11e-16 | PASS |
| `collinearity_tidy_pairs.csv` | 228 | 4.98e-13 | PASS |
| `collinearity_threshold_summary.csv` | 8 | 4.65e-13 | PASS |
| `collinearity_strongest_active7_pairs.csv` | 20 | 4.78e-13 | PASS |
| `collinearity_vif_active7.csv` | 14 | 4.23e-11 | PASS |

Predicted classes and fold memberships were checked for both years. SVM probability estimates are auxiliary outputs; they were not used to classify or choose parameters and are excluded from strict parity checking. No saved probability value was substituted for a newly calculated one.

## Fixed archived calls and RF refits

|   year | model   |   TN |   FP |   FN |   TP | matches_original_counts   |
|-------:|:--------|-----:|-----:|-----:|-----:|:--------------------------|
|   2005 | RF      |  251 |    4 |    7 |   37 | True                      |
|   2005 | SVM     |  248 |    7 |    6 |   38 | True                      |
|   2009 | RF      |  253 |    2 |    4 |   40 | True                      |
|   2009 | SVM     |  253 |    2 |    5 |   39 | True                      |

Both fixed archived SVM calls match the retained confusion counts. RF refit results are shown separately; even a count match is not proof of the complete historical fitted-object, tuning, or mapping provenance. The original manuscript metric table is always calculated from its explicitly labeled retained counts.

## Diagnostic equivalence

The collinearity and reference-label Moran diagnostics were originally run in Python. R equivalents reproduce the reported correlations, threshold counts, VIF values and reference-label Moran statistic within the stated tolerances. Each diagnostic uses the same 722 training IDs. R and NumPy permutation streams are not bit-identical; both recorded label-Moran p-values are 0.0001. Neighbor ties are resolved deterministically by row index in the R diagnostic.

Archived Python CSVs round many floats to 12 significant digits. Scaling constants can differ at that saved rounding precision; for example, the archived DEM mean can differ by about 4.18e-9 from the unrounded R value. Independent full-precision recalculation agrees within approximately 1e-12. No number was overwritten to hide rounding differences.

## Environment

See `outputs/sessionInfo.txt` for the complete R and package versions. The test environment uses R 4.3.3, randomForest 4.7-1.1 and e1071 1.7-14. Reproducibility outside this environment depends on input row order, package versions, R RNG settings and machine numeric behavior.

## What this test does not establish

It does not reproduce missing GEE/LandTrendr preprocessing, original whole-landscape maps, historical area totals, original feature-selection provenance or the historical Figure 7 importance ranking. It is not an independent validation of disturbance-agent labels, a feature-removal experiment, or a guarantee of generalization to new regions.
