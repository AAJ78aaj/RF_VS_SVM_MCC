# Clear R code for the RF–SVM manuscript

This package provides one commented analysis file, **`Manuscript_Analysis.R`**, with numbered sections. It consolidates the available code into a portable workflow for the manuscript’s binary classification comparison and supplementary reviewer analyses. It does not modify any manuscript, reviewer response, input CSV or original map.

## 1. Run the code

Extract `MPBB_R_Code.zip`. Keep the `data` folder beside `Manuscript_Analysis.R`.

In RStudio, open the script and select **Session → Set Working Directory → To Source File Location**. Install the two packages once:

```r
install.packages(c("randomForest", "e1071"))
```

Then run:

```r
source("Manuscript_Analysis.R")
```

From a terminal, you can instead run:

```bash
Rscript Manuscript_Analysis.R
```

All newly calculated results are saved under `outputs/`. The main run includes historical count verification, documented fixed-call refits, the primary spatial geometry, topography, predictor correlations and reference-label Moran’s I. To reproduce **all** of Table S5, set `RUN_SENSITIVITY <- TRUE` near the top of the script. Alternatively, in RStudio:

```r
Sys.setenv(MPBB_RUN_SENSITIVITY = "true")
source("Manuscript_Analysis.R")
```

The full sensitivity run evaluates four geometries in total: primary 5,000/1,000 units; 2,500/500 units; shifted primary origin `(0, 2500)`; and 7,500/1,500 units. These are **source-coordinate units**, not verified metres. Running the full analysis takes longer than the primary run because each additional geometry refits and tunes both classifiers.

## 2. Inputs and class coding

| Input | Year | Observations | Use |
|---|---|---:|---|
| `Book2005.csv` | 2005 | 1,021 | Binary RF–SVM comparison and spatial validation |
| `Book2009.csv` | 2009 | 1,021 | Binary RF–SVM comparison and spatial validation |
| `MPR_all2009.csv` | 2009 only | 1,445 | Supplementary topographic associations |

`PresenceAb = 1` means **infested** and `PresenceAb = 0` means **non-infested**. The classifier response is converted to a factor. The classification predictors, in their documented order, are:

```r
features <- c("band4", "NBR", "B54", "TCG", "DEM", "AnomB54", "DiffB54")
```

Coordinates are used to build spatial folds and buffers, not as RF/SVM predictors. `Percent_Cov` regression is not part of the binary comparison. The MPR logistic model is a separate supplementary analysis, not a replacement for classifier validation.

The seed-123 historical allocation produces **722 training** and **299 test** observations. It assigns a partition label to each row using probabilities 0.7/0.3; it is not a fixed-count split. Preserve the original CSV row order. The separately reported historical 714/307 diagnostic is recorded in the dataset-accounting output, but its memberships are not reconstructed or assumed to explain eight exclusions.

## 3. Read the script by section

| Section | Contents | Connection to the manuscript |
|---|---|---|
| 0 | Paths, packages, predictor list and spatial settings | Reproducible setup |
| 1 | Metrics, spatial blocks, buffering, fitting, nested tuning, block bootstrap and Moran helper functions | Methods 2.3–2.6 |
| 2 | Original confusion counts, corrected metrics, seed-123 partition, fixed archived SVM calls and a documented RF refit | Original Tables 2–4; supplementary Table S1 |
| 3 | Primary matched nested spatial validation in both years | Tables S3–S4 |
| 4 | 2009 aspect sectors, exact binomial intervals, chi-square, standardized logistic coefficients and residual Moran’s I | Table S6; topographic results |
| 5 | Training-only predictor correlations, auxiliary VIF diagnostics, reference-label Moran’s I and file-level class crosswalk | Table S2, Table S8 and spatial/reference reporting |
| 6 | Optional block-size/grid-origin sensitivity, input checksums and R session information | Table S5 and reproducibility records |

Table S7 is a narrative prior-work/reuse crosswalk, not a statistical calculation. Its content should remain in the manuscript; no numerical analysis creates a new originality claim.

## 4. Original results versus newly calculated results

### Original holdout summary

The original confusion counts are entered explicitly and labeled as **author-reported historical counts**:

| Year | Classifier | TN | FP | FN | TP |
|---|---|---:|---:|---:|---:|
| 2005 | RF | 251 | 4 | 7 | 37 |
| 2005 | SVM | 248 | 7 | 6 | 38 |
| 2009 | RF | 253 | 2 | 4 | 40 |
| 2009 | SVM | 253 | 2 | 5 | 39 |

`corrected_original_holdout_metrics.csv` recalculates OA, Kappa, MCC, precision, recall, specificity and F1 from these counts. This is a legitimate count-based verification, not a claim that the counts were regenerated from a newly recovered historical fitted object. In particular, SVM-2009 Kappa is approximately **0.904**, not 0.97.

### Documented fixed-call refits

The script separately fits the archived settings to the documented training partition:

```r
# RF, both years
randomForest(PresenceAb ~ ., data = training,
             ntree = 500, mtry = 4,
             importance = TRUE, proximity = TRUE)

# SVM, 2005
svm(PresenceAb ~ ., data = training,
    type = "C-classification", kernel = "radial",
    cost = 32, gamma = 0.1, scale = TRUE)

# SVM, 2009
svm(PresenceAb ~ ., data = training,
    type = "C-classification", kernel = "radial",
    cost = 8, gamma = 0.1, scale = TRUE)
```

`fixed_call_refit_metrics.csv` includes `matches_original_counts` for each model/year. It never substitutes a new RF fit’s counts for the original manuscript results. The archived RF call uses 500 trees; the historical narrative also mentions 600 trees. The package does not invent the fitted object or tuning record needed to resolve that distinction. A count match, if obtained, still does not establish the complete historical tuning provenance.

**Verified test result:** In the recorded R environment, the cleaned fixed calls reproduced all four original confusion matrices: RF and SVM in both 2005 and 2009. This is stronger than count-only metric verification, while remaining distinct from recovery of every historical tuning decision, original fitted object or raster-mapping workflow. See `VALIDATION_REPORT.md` for the exact counts and comparisons.

The archived scripts assessed classes using `caret::confusionMatrix()`. The cleaned script uses `table()` plus an explicit `metric_row()` function, so no caret dependency is required. These calculate metrics with positive class 1. MCC is not presented as a default caret output. The original scripts passed `epsilon = 0` to SVM; that does not define a C-classification loss and is omitted from the cleaned binary calls.

### Added spatial comparison

RF and SVM use the same seven predictors and outer test locations. The primary evaluation uses five outer spatial folds, a 5,000-unit square grid and a 1,000-unit training exclusion buffer. Three buffered inner folds select parameters by mean MCC within the outer training data:

- RF: 500 trees; `mtry = c(2, 4, 7)`.
- SVM: radial C-classification; `cost = c(0.25, 1, 4, 16)` and `gamma = c(1/14, 1/7, 2/7)`; scaling enabled; no class weights.

The 9,999-replicate paired block bootstrap resamples complete outer blocks. It conditions on existing out-of-fold predictions and does **not** refit models for each bootstrap replicate. One deterministic fold allocation is used per geometry. This is a validation-sensitivity analysis, not a spatial-autoregressive classifier or an isolated estimate of spatial bias.

### Topography and collinearity

The 2009 MPR model uses separately standardized DEM, slope, sine(aspect) and cosine(aspect). Coefficient/odds-ratio Wald intervals and exact binomial sector intervals are nominal, not spatially corrected. Residual Moran’s I documents remaining spatial dependence. ArcGIS extraction and binary coding do not, by themselves, establish independently observed outcome labels.

Predictor correlations use only the documented **722 training records per year**. Neither holdout values nor the response is used to calculate predictor correlations/VIF. Absolute-correlation thresholds are strict `> 0.6`, `> 0.7` and `> 0.8`. No predictor is removed or retrained on the basis of this diagnostic, and no performance benefit from screening is claimed. The active-seven and wider diagnostic-nine sets are distinguished; the missing NBR-anomaly column in the Book files is not invented or replaced by `AnomB54`.

## 5. Main output files

| Output | Meaning |
|---|---|
| `Table_S1_dataset_roles.csv` | Dataset roles, sizes and unresolved historical diagnostic membership |
| `archived_seed123_split_counts.csv` | 722/299 split and class counts |
| `archived_partition_locations.csv` | Point IDs, source coordinates, labels and partition membership |
| `corrected_original_holdout_metrics.csv` | Corrected metrics from the four retained original matrices |
| `fixed_call_refit_metrics.csv` | Separately calculated fixed-call refits and count-match flags |
| `fixed_call_holdout_predictions.csv` | Newly generated fixed-call predictions, not reconstructed original RF predictions |
| `spatial_oof_predictions.csv` | Paired RF/SVM outer-test predictions and fold/block IDs |
| `spatial_fold_accounting.csv` | Test/training class counts and buffer exclusions |
| `nested_tuning_results.csv` | All inner-grid scores, by year/model/outer fold |
| `spatial_performance_metrics.csv` | Primary spatial confusion counts and complete metrics |
| `spatial_block_bootstrap.csv` | Primary paired block-bootstrap estimates and intervals |
| `spatial_geometry_sensitivity.csv` | Full Table S5 results, created only when sensitivity is enabled |
| `mpr_aspect_rates_with_exact_ci.csv` | Sector counts, infestation proportions and exact intervals |
| `mpr_logistic_regression.csv` | Complete coefficient, SE, Wald CI and odds-ratio table |
| `mpr_predictor_scaling.csv` | Means and SDs used in the topographic model |
| `mpr_book_coordinate_overlap.csv` | Exact-coordinate overlap accounting |
| `reference_label_moran.csv` | Reference-label clustering, not a model predictor |
| `collinearity_tidy_pairs.csv` | Training-only Pearson/Spearman pair records |
| `collinearity_threshold_summary.csv` | Absolute-correlation threshold counts |
| `collinearity_vif_active7.csv` | Auxiliary linear redundancy diagnostic; not evidence of RF/SVM performance loss |
| `revised_analysis_summary.txt` | Spatial and topographic numerical summary |
| `sessionInfo.txt` | R and package versions |

Additional matrices and audit outputs are described in the code. CSV tables are editable data outputs; they do not automatically insert formatted tables into Word.

## 6. Provenance and limits

The `archive/` folder preserves the unmodified RF and SVM repository scripts, including their original paths, regression blocks and historical code issues. They are documentary sources, not automatically executed. It also preserves `revised_analysis.R`, the actual previously used script for the added spatial/topographic analysis. The clean main file carries its core algorithms and seeds forward with portable paths, comments, validation checks and additional documented exports.

The historical-metric, reference-label Moran and collinearity diagnostics were previously calculated in Python. This package supplies equivalent R calculations; it does not falsely describe those original diagnostics as having been run in R. R and NumPy random-number generators differ, so permutation-by-permutation correspondence is not expected. Discrete nearest-neighbor ties can also be resolved differently by software; the checked values and any differences are recorded in `VALIDATION_REPORT.md`.

This package starts from the **already extracted predictor CSVs**. It does not rerun GEE/LandTrendr, regenerate original full-landscape RF/SVM maps, calculate the historical infestation-area totals from rasters, or reproduce the historical Figure 7 importance ranking. The necessary original rasters, masks, fitted-object provenance and preprocessing records are not all supplied. No illustrative maps, new labels, independent infestation polygons or fabricated data are created.

Retain `sessionInfo.txt` and input checksums with any manuscript submission. Minor stochastic/software differences should be investigated, not silently replaced with the expected manuscript values.
