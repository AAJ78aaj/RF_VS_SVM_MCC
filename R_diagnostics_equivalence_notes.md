# Reference diagnostics: standalone base-R module

## Integration and scope

`diagnostics_module.R` defines exactly three top-level functions: `diagnose_knn`, `diagnose_moran`, and `run_reference_diagnostics`. It can be sourced or appended to the parent R file without renaming the parent's existing `knn_indices` or `moran_i_knn` helpers. It has no package dependencies, hard-coded input/output paths, global input variables, or execution on sourcing. Base/recommended R functions are used; no packages are attached. The only global state touched during execution is R's RNG machinery, whose kind and seed are restored on exit.

```r
source("diagnostics_module.R")
b05 <- read.csv("Book2005.csv", check.names = FALSE)
b09 <- read.csv("Book2009.csv", check.names = FALSE)
set.seed(123)
split_id <- sample(2, nrow(b05), replace = TRUE, prob = c(.7, .3))
diagnostics <- run_reference_diagnostics(b05, b09, split_id, out_dir)
```

The two Book arguments accept data.frames or individual CSV paths. `split_id` is a numeric vector in **Book2005 row order**, with `1` meaning training and `2` holdout. Training OIDs are matched separately to each year, so Book2009 need not share Book2005 row order. The function validates unique/nonmissing IDs, matching OID sets, correct years, binary labels, and matching locations/labels across years. Existing archived membership is used for independent validation, not as a hard-coded runtime dependency. The seed-123 allocation reproduces all archived training OIDs exactly.

These are diagnostics only. They never select/exclude predictors, refit RF/SVM, create lag predictors, correct predictions using Moran's I, or run a new suite of five significance tests. The auxiliary regressions used for VIF assess linear redundancy only. VIF is supplied as a secondary optional-use diagnostic; its inclusion here does **not** claim it was reported in the manuscript or used for original selection. Correlation/VIF do not establish classifier failure or reduced-feature performance.

## Training-only collinearity

The ordered sets are:

- `active7`: band4, NBR, B54, TCG, DEM, AnomB54, DiffB54.
- `diagnostic9`: band4, band5, NBR, B54, TCG, MeanB54, DEM, AnomB54, DiffB54.

All correlation, scaling, VIF, and column-audit calculations use only the same 722 archived training OIDs. Non-numeric/non-finite predictor values become missing for calculations, with their counts retained in the audit. Correlations use finite pairwise complete observations without imputation. Spearman uses pair-specific average ranks for ties. Undefined pairs retain statuses. Unique off-diagonal pairs are counted once, with strict unrounded thresholds `|r| > .6`, `>.7`, `>.8`.

The strongest-active7 file has ten distinct pairs per year ranked by `max(abs(Pearson), abs(Spearman))`, then predictor names in deterministic lexicographic order, as in Python. VIF uses complete-case active7 data, training means, sample SDs (`n-1`), an intercept and the other six predictors; `VIF = 1/(1-R²)`. Both supplied training designs have rank eight. No training rows or variables were removed for these supplied complete, finite data.

The CSVs provide only `NBR`, not an identifiable NBR-anomaly variable. No NBR anomaly is invented and `AnomB54` is not relabeled as NBR anomaly.

## Reference labels, spatial weights, and sample accounting

`reference_label_moran.csv` uses **all 1,021 reference labels separately for each year**, not training-only data, model residuals, or predictions. Euclidean distances use the original `Long`/`Lat` numeric **source projected-coordinate units**, not geographic degrees. No CRS, distance unit conversion, reprojection, or verified EPSG code is invented. Weights are directed eight-nearest-neighbor, row-standardized (`1/8` per neighbor); `S0=N`, so the Moran formula's `N/S0` multiplier is one. Self-neighbors are explicitly excluded. Equal distances are broken by source row index. The Python source used `scipy.spatial.cKDTree`; its tie selection can differ, potentially causing small I differences in other inputs. No substantive I difference arose here.

Each year resets `set.seed(20261104)` using Mersenne-Twister/Inversion/Rejection and runs 9,999 label permutations with fixed neighbors. The two-sided Monte Carlo p is `(1 + count(abs(I_perm) >= abs(I_obs))) / 10000`. This R RNG stream is **not bit-identical to NumPy's** `default_rng`; exact permutation sequences are not claimed equivalent. Here both implementations return `p = 0.0001`, the minimum attainable p, rather than `p < 0.0001`. Constant labels give undefined I/p. Equal annual results reflect identical locations and labels, not independent annual reobservation.

`reference_sample_accounting.csv` records actual reference/train/holdout class counts and the nominal Cochran calculation `qnorm(.975)^2 * p * (1-p) / .05^2`, without a finite-population correction. It gives `205.122295293509`, rounded upward to **206**, not 205. This nominal independent-sampling prevalence calculation does not establish recall precision, classifier-comparison power, spatial independence, or representative coverage. It does not reconstruct or confuse the historical 714/307 diagnostic with the reproducible 722/299 allocation.

`reference_code_crosswalk.csv` is record-level and preserves each year's OID and raw Classcode, Classname, Classvalue, and PresenceAb. `Classvar` is inspected by exact column name: it is absent in both supplied CSVs, so the output contains `NA` and `Classvar_present=FALSE`. Absence codes are 10/10 and binary label 0; presence codes are 20/20 and label 1. The source spelling variants `Presence_Dead` and `Presence Dead` remain unchanged; no synthetic identifier, class recoding, or `Classvar` alias is created.

## Ordered output schemas

All CSVs have no row-number column. Names and column order of the original collinearity and Moran outputs are preserved. CSV quoting, logical capitalization, numeric display precision, and scientific notation can differ from pandas but not the data/schema. Numeric results are not rounded before thresholding or ranking.

| Filename | Columns, in order |
| --- | --- |
| `correlation_{2005,2009}_{active7,diagnostic9}_{pearson,spearman}.csv` (eight files) | `predictor`, followed by the predictors in the corresponding ordered set above |
| `collinearity_tidy_pairs.csv` | `year,predictor_set,method,predictor_1,predictor_2,n_pair,correlation,abs_correlation,status` |
| `collinearity_threshold_summary.csv` | `year,predictor_set,method,training_n,total_unique_pairs,valid_unique_pairs,undefined_pairs,abs_r_gt_0_6,abs_r_gt_0_7,abs_r_gt_0_8,strongest_pair,strongest_correlation` |
| `collinearity_strongest_active7_pairs.csv` | `year,rank,predictor_1,predictor_2,n_pair,pearson,spearman,max_abs_correlation` |
| `collinearity_vif_active7.csv` | `year,predictor,n_complete_case,R2_auxiliary,VIF,design_rank,status` |
| `collinearity_vif_training_scaling.csv` | `year,predictor,n_complete_case,training_mean,training_sample_sd` |
| `collinearity_column_audit.csv` | `year,predictor,training_n,raw_missing_n,non_numeric_n,nonfinite_nonmissing_n,finite_n,finite_unique_n,constant` |
| `reference_label_moran.csv` | `year,N,k,weights,Moran_I,two_sided_permutation_p,permutations,seed` |
| `reference_sample_accounting.csv` | `year,N,positive_n,negative_n,training_n,training_positive_n,training_negative_n,holdout_n,holdout_positive_n,holdout_negative_n,prevalence,confidence_level,margin_of_error,z_critical,cochran_n_unrounded,cochran_n_ceiling` |
| `reference_code_crosswalk.csv` | `year,OID_,Classcode,Classname,Classvalue,Classvar,PresenceAb,partition,Classvar_present` |
| `reference_diagnostic_membership.csv` | `year,OID_,partition` |

Tidy rows follow the Python sort: year, predictor set, method, descending absolute correlation, with stable combination order for ties. Summaries follow year, active7/diagnostic9, Pearson/Spearman. VIF/scaling follow year and active7 order. The function invisibly returns the ten tabular outputs as a named list; matrices are saved separately. It does not generate the old JSON/prose audit reports or corrected model-performance tables.

## Validation completed

Isolated tests and generated outputs are under `/home/ubuntu/revision_work/r_code_package/diagnostics_test_outputs/`. The completed log is `validation_log.txt`; the reusable R test is `validate_diagnostics.R`. Run it with:

```sh
Rscript /home/ubuntu/revision_work/r_code_package/diagnostics_test_outputs/validate_diagnostics.R
```

All tests passed. They check ordered schemas and row order, archived training OIDs, all eight correlation matrices, tidy pairs, threshold summaries, top-ten ranking, VIF, scaling, audit counts, Moran, accounting, record codes, RNG-state preservation, deterministic distance ties, duplicate-coordinate self exclusion, and k=1 helper dimensions. They also perturb every holdout predictor and reverse Book2009 row order: training correlation/VIF/scaling/audit results remain unchanged. Protected inputs and original Python scripts have unchanged hashes.

| Check | Maximum absolute difference / finding |
| --- | --- |
| Pearson/Spearman matrices and tidy pairs vs archived Python CSVs | `4.97990537695614e-13`, below `1e-10` |
| VIF vs archived Python CSV | `4.23021617734776e-11`, below `1e-10` |
| Full-precision scaling vs independently recomputed pandas means/sample SDs | `9.09494701772928e-13`, below `1e-10` |
| Moran I, R minus archived Python | `-1.11022302462516e-16` in both years |
| Moran I, R minus rounded manuscript `0.836` | `-0.000202740769485010` in both years; both round to `0.836` |
| Actual sample accounting, each year | N=1,021: 162 positive/859 negative; training=722: 118/604; holdout=299: 44/255 |

**Archived scaling precision caveat:** Python's `save_csv` used `float_format='%.12g'`. The saved DEM mean therefore differs from the unrounded R mean by `4.18003764934838e-9` and its SD by about `3.19e-10`. This is file-display rounding, not a computation discrepancy. Direct absolute `1e-10` comparison to that rounded scaling file cannot pass honestly. At the archived 12-significant-digit precision the scaling matches; independently recomputed full-precision Python scaling passes strict `1e-10`. The module retains full precision rather than distorting its calculations to match archival display. `validation_numeric_comparison.csv` distinguishes the checks explicitly.

Threshold counts match Python exactly:

| Year | Set | Method | Unique pairs | >.6 | >.7 | >.8 |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| 2005 | active7 | Pearson | 21 | 3 | 3 | 2 |
| 2005 | active7 | Spearman | 21 | 3 | 2 | 2 |
| 2005 | diagnostic9 | Pearson | 36 | 8 | 8 | 6 |
| 2005 | diagnostic9 | Spearman | 36 | 8 | 7 | 5 |
| 2009 | active7 | Pearson | 21 | 4 | 2 | 2 |
| 2009 | active7 | Spearman | 21 | 3 | 2 | 1 |
| 2009 | diagnostic9 | Pearson | 36 | 9 | 7 | 7 |
| 2009 | diagnostic9 | Spearman | 36 | 8 | 7 | 2 |

The strongest pair is NBR/B54 throughout: 2005 Pearson `-0.978226947173`, Spearman `-0.964782754058`; 2009 Pearson `-0.973485063525`, Spearman `-0.937064858958`. Both years' top-three active7 pairs are NBR/B54, band4/TCG, and AnomB54/DiffB54. Full coefficients and all ten distinct pairs are retained in the CSVs and validation log. Appreciable redundancy is confirmed; a universal `|r| <= .6` statement about the retained active7 set is not supported by these diagnostics.
