# ============================================================================
# RF versus SVM for mountain pine beetle: manuscript analysis code
# Read README.md first. This is a documented, portable consolidation, not the
# original historical author script. Original scripts are retained in archive/.
#
# What is reproducible here:
#   A. Metrics calculated from the FOUR ORIGINAL confusion matrices.
#   B. Fixed archived SVM binary calls and a documented RF refit for comparison.
#   C. Matched nested buffered spatial validation and paired uncertainty.
#   D. The separate 2009 MPR topographic analysis.
#   E. Training-only collinearity and reference-label spatial diagnostics.
#
# IMPORTANT: Recomputing original RF metrics from published counts is not proof
# of the exact historical RF fit, its tuning, or raster-map reproduction.
# Collinearity and label-Moran diagnostics were originally calculated in Python;
# their R equivalents are supplied below and numerical agreement is checked.
#
# RUN IN RSTUDIO: Set the working directory to the extracted MPBB_R_Code folder,
# then open this file and click Source. Or use the terminal:
#   Rscript Manuscript_Analysis.R
# Install once, if needed: install.packages(c("randomForest", "e1071"))
# ============================================================================

# 0. CONFIGURATION ------------------------------------------------------------
# Paths work on Windows, macOS, and Linux. No original input is modified.
# When using RStudio, setwd("C:/your/path/MPBB_R_Code") before sourcing this file.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
PROJECT_DIR <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE))
} else getwd()
INPUT_DIR <- file.path(PROJECT_DIR, "data")
OUTPUT_DIR <- Sys.getenv("MPBB_OUTPUT_DIR", unset = file.path(PROJECT_DIR, "outputs"))
RUN_SENSITIVITY <- tolower(Sys.getenv("MPBB_RUN_SENSITIVITY", "false")) %in% c("true", "1", "yes")
# RUN_SENSITIVITY <- TRUE  # Uncomment to reproduce all four Table S5 geometries.

needed <- c("randomForest", "e1071")
missing_packages <- needed[!vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) stop("Install missing packages first: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({library(randomForest); library(e1071)})
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
input_dir <- INPUT_DIR
out_dir <- OUTPUT_DIR
features <- c("band4", "NBR", "B54", "TCG", "DEM", "AnomB54", "DiffB54")
# PresenceAb: 1 = infested; 0 = non-infested. This is CLASSIFICATION, not
# Percent_Cov regression. Coordinates are NOT included as model predictors.
# Long/Lat are source-coordinate fields, not assumed geographic degrees/metres.
block_size <- 5000
buffer_size <- 1000
outer_shift_x <- 0
outer_shift_y <- 0
outer_k <- 5
inner_k <- 3
rf_mtry_grid <- c(2, 4, 7)
svm_cost_grid <- c(0.25, 1, 4, 16)
svm_gamma_grid <- c(1/14, 1/7, 2/7)

# 1. REUSABLE FUNCTIONS --------------------------------------------------------
# Confusion-count convention: TP/FP refer to predictions of positive class 1.
# Binary confusion-matrix exports below use reference as rows and prediction
# as columns. MCC uses all four counts; Kappa uses the observed class marginals.
# The MCC denominator uses double precision to avoid integer overflow.
metric_row <- function(obs, pred) {
  obs <- factor(as.character(obs), levels = c("0", "1"))
  pred <- factor(as.character(pred), levels = c("0", "1"))
  tn <- sum(obs == "0" & pred == "0")
  fp <- sum(obs == "0" & pred == "1")
  fn <- sum(obs == "1" & pred == "0")
  tp <- sum(obs == "1" & pred == "1")
  n <- tn + fp + fn + tp
  oa <- (tp + tn) / n
  precision <- if ((tp + fp) == 0) NA else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0) NA else tn / (tn + fp)
  f1 <- if (is.na(precision) || is.na(recall) || precision + recall == 0) NA else 2 * precision * recall / (precision + recall)
  pe <- (((tp + fp) * (tp + fn)) + ((fn + tn) * (fp + tn))) / n^2
  kappa <- if (abs(1 - pe) < .Machine$double.eps) NA else (oa - pe) / (1 - pe)
  den <- sqrt(as.numeric(tp + fp) * as.numeric(tp + fn) * as.numeric(tn + fp) * as.numeric(tn + fn))
  mcc <- if (den == 0) NA else (tp * tn - fp * fn) / den
  data.frame(TN = tn, FP = fp, FN = fn, TP = tp, N = n, OA = oa, Kappa = kappa,
             MCC = mcc, Precision = precision, Recall = recall, Specificity = specificity, F1 = f1)
}

make_blocks <- function(dat, size, shift_x = 0, shift_y = 0) {
  bx <- floor((dat$Long - min(dat$Long) + shift_x) / size)
  by <- floor((dat$Lat - min(dat$Lat) + shift_y) / size)
  paste(bx, by, sep = "_")
}

assign_folds <- function(block, y, k) {
  tab <- aggregate(data.frame(n = rep(1, length(block)), pos = as.integer(y == 1)),
                   by = list(block = block), FUN = sum)
  tab <- tab[order(-tab$pos, -tab$n, tab$block), ]
  fold_n <- rep(0, k); fold_pos <- rep(0, k)
  tab$fold <- NA_integer_
  for (i in seq_len(nrow(tab))) {
    score <- (fold_pos / max(1, sum(tab$pos))) + (fold_n / max(1, sum(tab$n)))
    target <- which(score == min(score))
    # break ties deterministically by smallest fold index
    f <- target[1]
    tab$fold[i] <- f
    fold_n[f] <- fold_n[f] + tab$n[i]
    fold_pos[f] <- fold_pos[f] + tab$pos[i]
  }
  tab$fold[match(block, tab$block)]
}

buffer_train_idx <- function(dat, test_idx, candidate_idx, buffer) {
  if (buffer <= 0) return(candidate_idx)
  te <- as.matrix(dat[test_idx, c("Long", "Lat")])
  tr <- as.matrix(dat[candidate_idx, c("Long", "Lat")])
  d2 <- outer(tr[, 1], te[, 1], "-")^2 + outer(tr[, 2], te[, 2], "-")^2
  candidate_idx[apply(d2, 1, min) > buffer^2]
}

fit_predict_rf <- function(train, test, mtry, seed) {
  set.seed(seed)
  fit <- randomForest(x = train[, features, drop = FALSE], y = factor(train$PresenceAb, levels = c(0, 1)),
                      ntree = 500, mtry = mtry, importance = TRUE)
  pred <- predict(fit, test[, features, drop = FALSE], type = "response")
  prob <- predict(fit, test[, features, drop = FALSE], type = "prob")[, "1"]
  list(pred = pred, prob = prob, fit = fit)
}

fit_predict_svm <- function(train, test, cost, gamma) {
  fit <- svm(x = train[, features, drop = FALSE], y = factor(train$PresenceAb, levels = c(0, 1)),
             type = "C-classification", kernel = "radial", cost = cost, gamma = gamma,
             scale = TRUE, class.weights = NULL, probability = TRUE)
  pred <- predict(fit, test[, features, drop = FALSE], probability = TRUE)
  prob <- attr(pred, "probabilities")[, "1"]
  list(pred = pred, prob = prob, fit = fit)
}

inner_score <- function(train, model, params, year_index, outer_fold) {
  train$block_inner <- make_blocks(train, block_size, shift_x = block_size / 2, shift_y = block_size / 2)
  k_use <- min(inner_k, length(unique(train$block_inner)))
  if (k_use < 2) return(NA_real_)
  train$fold_inner <- assign_folds(train$block_inner, train$PresenceAb, k_use)
  scores <- numeric(0)
  for (f in seq_len(k_use)) {
    i_test <- which(train$fold_inner == f)
    candidates <- which(train$fold_inner != f)
    i_train <- buffer_train_idx(train, i_test, candidates, buffer_size)
    if (length(i_train) < 30 || length(unique(train$PresenceAb[i_train])) < 2 || length(unique(train$PresenceAb[i_test])) < 2) next
    ans <- tryCatch({
      if (model == "RF") fit_predict_rf(train[i_train, ], train[i_test, ], params$mtry,
                                          seed = 100000 + year_index * 1000 + outer_fold * 100 + f)
      else fit_predict_svm(train[i_train, ], train[i_test, ], params$cost, params$gamma)
    }, error = function(e) NULL)
    if (!is.null(ans)) scores <- c(scores, metric_row(train$PresenceAb[i_test], ans$pred)$MCC)
  }
  if (!length(scores)) NA_real_ else mean(scores, na.rm = TRUE)
}

tune_model <- function(train, model, year_index, outer_fold) {
  if (model == "RF") {
    grid <- data.frame(mtry = rf_mtry_grid)
    grid$inner_MCC <- vapply(seq_len(nrow(grid)), function(i) {
      inner_score(train, "RF", grid[i, , drop = FALSE], year_index, outer_fold)
    }, numeric(1))
    if (all(is.na(grid$inner_MCC))) stop("No valid inner-fold tuning scores.")
    best <- grid[which.max(replace(grid$inner_MCC, is.na(grid$inner_MCC), -Inf)), , drop = FALSE]
  } else {
    grid <- expand.grid(cost = svm_cost_grid, gamma = svm_gamma_grid)
    grid$inner_MCC <- vapply(seq_len(nrow(grid)), function(i) {
      inner_score(train, "SVM", grid[i, , drop = FALSE], year_index, outer_fold)
    }, numeric(1))
    if (all(is.na(grid$inner_MCC))) stop("No valid inner-fold tuning scores.")
    best <- grid[which.max(replace(grid$inner_MCC, is.na(grid$inner_MCC), -Inf)), , drop = FALSE]
  }
  list(grid = grid, best = best)
}

run_spatial_nested <- function(dat, year_label, year_index) {
  dat$block_outer <- make_blocks(dat, block_size, shift_x = outer_shift_x, shift_y = outer_shift_y)
  dat$fold_outer <- assign_folds(dat$block_outer, dat$PresenceAb, outer_k)
  oof <- list(); tuning <- list(); fold_info <- list()
  for (f in seq_len(outer_k)) {
    message(sprintf("Spatial validation: year %s, outer fold %d/%d", year_label, f, outer_k))
    i_test <- which(dat$fold_outer == f)
    candidates <- which(dat$fold_outer != f)
    i_train <- buffer_train_idx(dat, i_test, candidates, buffer_size)
    if (length(i_train) < 30 || length(unique(dat$PresenceAb[i_train])) < 2 || length(unique(dat$PresenceAb[i_test])) < 2) {
      stop(sprintf("Invalid outer fold %s for %s after buffer: train=%s test=%s", f, year_label, length(i_train), length(i_test)))
    }
    tr <- dat[i_train, , drop = FALSE]
    te <- dat[i_test, , drop = FALSE]
    trf <- tune_model(tr, "RF", year_index, f)
    tsvm <- tune_model(tr, "SVM", year_index, f)
    rf <- fit_predict_rf(tr, te, trf$best$mtry[1], seed = 10000 + year_index * 100 + f)
    sv <- fit_predict_svm(tr, te, tsvm$best$cost[1], tsvm$best$gamma[1])
    oof[[f]] <- data.frame(year = year_label, OID = dat$OID_[i_test], Long = dat$Long[i_test], Lat = dat$Lat[i_test],
                           observed = dat$PresenceAb[i_test], RF_prediction = as.integer(as.character(rf$pred)), RF_probability = rf$prob,
                           SVM_prediction = as.integer(as.character(sv$pred)), SVM_probability = sv$prob,
                           outer_block = dat$block_outer[i_test], outer_fold = f)
    tuning[[f]] <- rbind(
      data.frame(year = year_label, outer_fold = f, model = "RF", mtry = trf$grid$mtry,
                 cost = NA_real_, gamma = NA_real_, inner_MCC = trf$grid$inner_MCC),
      data.frame(year = year_label, outer_fold = f, model = "SVM", mtry = NA_real_,
                 cost = tsvm$grid$cost, gamma = tsvm$grid$gamma, inner_MCC = tsvm$grid$inner_MCC)
    )
    fold_info[[f]] <- data.frame(year = year_label, outer_fold = f,
                                 test_n = length(i_test), test_negative = sum(te$PresenceAb == 0), test_positive = sum(te$PresenceAb == 1),
                                 train_n_after_buffer = length(i_train), train_negative = sum(tr$PresenceAb == 0), train_positive = sum(tr$PresenceAb == 1),
                                 excluded_by_buffer = length(candidates) - length(i_train),
                                 test_blocks = length(unique(te$block_outer)))
  }
  list(oof = do.call(rbind, oof), tuning = do.call(rbind, tuning), folds = do.call(rbind, fold_info), data = dat)
}

block_bootstrap <- function(pred, blocks, B = 9999, seed = 451) {
  set.seed(seed)
  ub <- unique(blocks)
  get_stats <- function(rf, sv) {
    c(RF_MCC = rf$MCC, RF_Precision = rf$Precision, RF_Recall = rf$Recall,
      RF_Specificity = rf$Specificity, RF_F1 = rf$F1,
      SVM_MCC = sv$MCC, SVM_Precision = sv$Precision, SVM_Recall = sv$Recall,
      SVM_Specificity = sv$Specificity, SVM_F1 = sv$F1,
      Delta_MCC_RF_minus_SVM = rf$MCC - sv$MCC,
      Delta_Precision_RF_minus_SVM = rf$Precision - sv$Precision,
      Delta_Recall_RF_minus_SVM = rf$Recall - sv$Recall)
  }
  first_rf <- metric_row(pred$observed, pred$RF_prediction)
  first_sv <- metric_row(pred$observed, pred$SVM_prediction)
  estimate <- get_stats(first_rf, first_sv)
  out <- matrix(NA_real_, B, length(estimate))
  colnames(out) <- names(estimate)
  for (b in seq_len(B)) {
    samp <- sample(ub, length(ub), replace = TRUE)
    idx <- unlist(lapply(samp, function(z) which(blocks == z)), use.names = FALSE)
    rf <- metric_row(pred$observed[idx], pred$RF_prediction[idx])
    sv <- metric_row(pred$observed[idx], pred$SVM_prediction[idx])
    out[b, ] <- get_stats(rf, sv)
  }
  data.frame(metric = names(estimate), estimate = as.numeric(estimate),
             CI_low = apply(out, 2, quantile, probs = 0.025, na.rm = TRUE),
             CI_high = apply(out, 2, quantile, probs = 0.975, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

knn_indices <- function(x, y, k = 8) {
  n <- length(x); d2 <- outer(x, x, "-")^2 + outer(y, y, "-")^2; diag(d2) <- Inf
  t(apply(d2, 1, function(z) order(z)[seq_len(k)]))
}
moran_i_knn <- function(z, neighbors) {
  zz <- z - mean(z)
  sum(zz * rowMeans(matrix(zz[neighbors], nrow = nrow(neighbors), ncol = ncol(neighbors)))) / sum(zz^2)
}

# 2. ORIGINAL HOLDOUT: COUNT VERIFICATION AND DOCUMENTED REFITS ----------------
# This is the seed-123 probabilistic allocation in the archived R scripts.
# replace=TRUE assigns one partition label to each row; it does not duplicate
# observations. A nominal 70/30 probability does not imply exact 714/307 counts.
# Preserve row order: sorting input changes which records receive each label.
b05 <- read.csv(file.path(input_dir, "Book2005.csv"), check.names = FALSE)
b09 <- read.csv(file.path(input_dir, "Book2009.csv"), check.names = FALSE)
mpr_check <- read.csv(file.path(input_dir, "MPR_all2009.csv"), check.names = FALSE)
check_book <- function(d, year) {
  required <- c("OID_", "year", "Long", "Lat", "PresenceAb", features)
  stopifnot(all(required %in% names(d)), nrow(d) == 1021,
            all(d$year == year), !anyDuplicated(d$OID_),
            all(d$PresenceAb %in% c(0, 1)))
  stopifnot(all(vapply(d[, c(features, "Long", "Lat")],
                       function(x) is.numeric(x) && all(is.finite(x)), logical(1))))
}
check_book(b05, 2005); check_book(b09, 2009)
stopifnot(identical(b05$OID_, b09$OID_), identical(b05$Long, b09$Long),
          identical(b05$Lat, b09$Lat), identical(b05$PresenceAb, b09$PresenceAb))
stopifnot(nrow(mpr_check) == 1445, all(mpr_check$Year == 2009),
          all(mpr_check$PresenceAb %in% c(0, 1)))
set.seed(123)
archived_split <- sample(2, nrow(b05), replace = TRUE, prob = c(0.7, 0.3))
archived_counts <- data.frame(
  partition = c("training", "test"),
  n = c(sum(archived_split == 1), sum(archived_split == 2)),
  negative = c(sum(b05$PresenceAb[archived_split == 1] == 0),
               sum(b05$PresenceAb[archived_split == 2] == 0)),
  positive = c(sum(b05$PresenceAb[archived_split == 1] == 1),
               sum(b05$PresenceAb[archived_split == 2] == 1)))
stopifnot(identical(archived_counts$n, c(722L, 299L)))
write.csv(archived_counts, file.path(out_dir, "archived_seed123_split_counts.csv"), row.names = FALSE)
locations <- data.frame(b05[, c("OID_", "Long", "Lat", "PresenceAb")],
                        partition = ifelse(archived_split == 1, "training", "test"))
write.csv(locations, file.path(out_dir, "archived_partition_locations.csv"), row.names = FALSE)

# These are the retained author-reported counts, NOT predictions reconstructed
# by this script. The manuscript's corrected original metrics use these counts.
original_counts <- data.frame(
  year = c(2005, 2005, 2009, 2009), model = c("RF", "SVM", "RF", "SVM"),
  TN = c(251, 248, 253, 253), FP = c(4, 7, 2, 2),
  FN = c(7, 6, 4, 5), TP = c(37, 38, 40, 39))
metrics_from_counts <- function(tn, fp, fn, tp) {
  observed <- c(rep(0, tn + fp), rep(1, fn + tp))
  predicted <- c(rep(0, tn), rep(1, fp), rep(0, fn), rep(1, tp))
  metric_row(observed, predicted)
}
historical_metrics <- do.call(rbind, lapply(seq_len(nrow(original_counts)), function(i) {
  x <- original_counts[i, ]
  data.frame(year = x$year, model = x$model,
             metrics_from_counts(x$TN, x$FP, x$FN, x$TP))
}))
# SVM 2009 Kappa = approximately 0.904, NOT 0.97.
write.csv(historical_metrics, file.path(out_dir, "corrected_original_holdout_metrics.csv"), row.names = FALSE)
print(historical_metrics)

# Clean fixed-call workflow using the documented settings. RF is labelled as
# a refit rather than guaranteed recovery of the author's historical RF object.
# Full historical feature-selection/tuning provenance is not reconstructed.
refit_metrics <- list(); refit_predictions <- list()
for (year in c(2005, 2009)) {
  d <- if (year == 2005) b05 else b09
  d$PresenceAb <- factor(d$PresenceAb, levels = c(0, 1))
  set.seed(123)
  split_id <- sample(2, nrow(d), replace = TRUE, prob = c(0.7, 0.3))
  training <- d[split_id == 1, c("PresenceAb", features)]
  testing <- d[split_id == 2, c("PresenceAb", features)]
  rf_fixed <- randomForest(PresenceAb ~ ., data = training, ntree = 500,
                          mtry = 4, importance = TRUE, proximity = TRUE)
  svm_fixed <- svm(PresenceAb ~ ., data = training,
                   type = "C-classification", kernel = "radial", scale = TRUE,
                   cost = if (year == 2005) 32 else 8, gamma = 0.1)
  # epsilon controls regression, not C-classification. It is not a tuned
  # classification tolerance here. The original scripts passed epsilon=0.
  rf_class <- predict(rf_fixed, testing)
  svm_class <- predict(svm_fixed, testing)
  for (model in c("RF", "SVM")) {
    pred <- if (model == "RF") rf_class else svm_class
    m <- metric_row(testing$PresenceAb, pred)
    expected <- original_counts[original_counts$year == year & original_counts$model == model, ]
    matches <- all(as.numeric(m[1, c("TN", "FP", "FN", "TP")]) ==
                   as.numeric(expected[1, c("TN", "FP", "FN", "TP")]))
    refit_metrics[[length(refit_metrics) + 1]] <- data.frame(year = year, model = model,
                                                          matches_original_counts = matches, m)
    cm <- table(reference = factor(testing$PresenceAb, levels = c(0, 1)),
                prediction = factor(pred, levels = c(0, 1)))
    write.csv(as.data.frame.matrix(cm), file.path(out_dir, paste0("refit_", model, "_", year, "_confusion_matrix.csv")))
  }
  refit_predictions[[as.character(year)]] <- data.frame(
    year = year, OID_ = d$OID_[split_id == 2], observed = testing$PresenceAb,
    RF_prediction = rf_class, SVM_prediction = svm_class)
  saveRDS(rf_fixed, file.path(out_dir, paste0("documented_RF_refit_", year, ".rds")))
  saveRDS(svm_fixed, file.path(out_dir, paste0("fixed_archived_SVM_call_", year, ".rds")))
}
write.csv(do.call(rbind, refit_metrics), file.path(out_dir, "fixed_call_refit_metrics.csv"), row.names = FALSE)
write.csv(do.call(rbind, refit_predictions), file.path(out_dir, "fixed_call_holdout_predictions.csv"), row.names = FALSE)
# Do not substitute new RF refit counts for the retained manuscript counts.

# Table S1 records dataset roles. The 714/307 diagnostic is only a documented
# historical description; its memberships cannot be recovered from these CSVs.
S1 <- data.frame(
  dataset = c("Book reference", "Seed-123 training", "Seed-123 test", "Historical spatial diagnostic",
              "New spatial out-of-fold comparison", "MPR_all2009"),
  years = c("2005 and 2009", "2005 and 2009", "2005 and 2009", "2009-designated",
            "2005 and 2009", "2009 only"),
  N = c("1021 each", "722 each", "299 each", "714 training / 307 diagnostic", "1021 each", "1445"),
  positive = c("162", "118", "44", "49 in diagnostic (historically reported)", "162", "120"),
  negative = c("859", "604", "255", "258 in diagnostic (historically reported)", "859", "1325"),
  role = c("Original comparison and added spatial validation; shared locations/labels",
           "Historical fixed-call model training", "Original reported holdout matrices",
           "No retained membership record; not assumed to be 307 minus 8",
           "New matched nested validation; every point tested once",
           "ArcGIS-extracted random points for supplementary topography, not independent classifier validation"))
write.csv(S1, file.path(out_dir, "Table_S1_dataset_roles.csv"), row.names = FALSE)

# 3. MATCHED NESTED BUFFERED SPATIAL VALIDATION ---------------------------------
# Same outer test locations and seven predictors for RF/SVM. Inner spatial CV
# tunes only on outer-training data. The bootstrap resamples complete blocks,
# conditioning on fitted predictions: it DOES NOT refit models per replicate.
# The design uses outcome counts for balancing blocks, but no test predictors
# are used to tune models. One deterministic allocation is used per geometry.
# New, fully documented matched nested spatial CV.
r05 <- run_spatial_nested(b05, "2005", 1)
r09 <- run_spatial_nested(b09, "2009", 2)
oof <- rbind(r05$oof, r09$oof)
tuning <- rbind(r05$tuning, r09$tuning)
folds <- rbind(r05$folds, r09$folds)
write.csv(oof, file.path(out_dir, "spatial_oof_predictions.csv"), row.names = FALSE)
write.csv(tuning, file.path(out_dir, "nested_tuning_results.csv"), row.names = FALSE)
write.csv(folds, file.path(out_dir, "spatial_fold_accounting.csv"), row.names = FALSE)

perf <- do.call(rbind, lapply(split(oof, oof$year), function(x) {
  rf <- metric_row(x$observed, x$RF_prediction); rf$model <- "RF"; rf$year <- unique(x$year)
  sv <- metric_row(x$observed, x$SVM_prediction); sv$model <- "SVM"; sv$year <- unique(x$year)
  rbind(rf, sv)
}))
perf <- perf[, c("year", "model", "N", "TN", "FP", "FN", "TP", "OA", "Kappa", "MCC", "Precision", "Recall", "Specificity", "F1")]
write.csv(perf, file.path(out_dir, "spatial_performance_metrics.csv"), row.names = FALSE)

boot <- do.call(rbind, lapply(split(oof, oof$year), function(x) {
  z <- block_bootstrap(x, x$outer_block, B = 9999, seed = 450 + as.integer(unique(x$year)))
  cbind(year = unique(x$year), z)
}))
write.csv(boot, file.path(out_dir, "spatial_block_bootstrap.csv"), row.names = FALSE)


# 4. SEPARATE 2009 TOPOGRAPHIC ANALYSIS -----------------------------------------
# MPR: random spatial points with attributes extracted in ArcGIS.
# This section is supplementary descriptive topography, not classifier testing.
# Exact binomial and Wald intervals are nominal: not spatially corrected.
# MPR topographic analysis, reported explicitly as a supplemental overlapping sample.
mpr <- read.csv(file.path(input_dir, "MPR_all2009.csv"), check.names = FALSE)
stopifnot(all(mpr$Year == 2009), nrow(mpr) == 1445,
          all(vapply(mpr[, c("DEM", "Slope", "Aspect", "Long", "Lat")],
                     function(x) is.numeric(x) && all(is.finite(x)), logical(1))),
          all(mpr$Aspect >= 0 & mpr$Aspect <= 360), all(mpr$Slope >= 0))
aspect_label <- function(a) {
  labs <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  labs[(floor(((a + 22.5) %% 360) / 45) + 1)]
}
mpr$aspect_class <- factor(aspect_label(mpr$Aspect), levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))
tab <- table(mpr$aspect_class, mpr$PresenceAb)
chisq_asymptotic <- chisq.test(tab, correct = FALSE)
set.seed(20261031)
chisq_montecarlo <- chisq.test(tab, correct = FALSE, simulate.p.value = TRUE, B = 99999)
aspect_summary <- data.frame(aspect = rownames(tab), healthy = tab[, "0"], infested = tab[, "1"], total = rowSums(tab))
aspect_summary$rate <- aspect_summary$infested / aspect_summary$total
ci <- t(mapply(function(x, n) {
  z <- binom.test(x, n)$conf.int
  c(z[1], z[2])
}, aspect_summary$infested, aspect_summary$total))
aspect_summary$CI_low <- ci[, 1]
aspect_summary$CI_high <- ci[, 2]
write.csv(aspect_summary, file.path(out_dir, "mpr_aspect_rates_with_exact_ci.csv"), row.names = FALSE)

mpr$z_DEM <- as.numeric(scale(mpr$DEM))
mpr$z_Slope <- as.numeric(scale(mpr$Slope))
mpr$z_aspect_sin <- as.numeric(scale(sin(mpr$Aspect * pi / 180)))
mpr$z_aspect_cos <- as.numeric(scale(cos(mpr$Aspect * pi / 180)))
aspect_scaling <- data.frame(variable = c("DEM", "Slope", "sin(aspect)", "cos(aspect)"),
                             mean = c(mean(mpr$DEM), mean(mpr$Slope), mean(sin(mpr$Aspect * pi / 180)), mean(cos(mpr$Aspect * pi / 180))),
                             sd = c(sd(mpr$DEM), sd(mpr$Slope), sd(sin(mpr$Aspect * pi / 180)), sd(cos(mpr$Aspect * pi / 180))))
write.csv(aspect_scaling, file.path(out_dir, "mpr_predictor_scaling.csv"), row.names = FALSE)
topo_fit <- glm(PresenceAb ~ z_DEM + z_Slope + z_aspect_sin + z_aspect_cos, family = binomial(), data = mpr)
co <- summary(topo_fit)$coefficients
ci_logit <- suppressMessages(confint.default(topo_fit))
topo_out <- data.frame(term = rownames(co), estimate = co[, "Estimate"], std_error = co[, "Std. Error"], z = co[, "z value"], p_value = co[, "Pr(>|z|)"],
                       CI_low = ci_logit[, 1], CI_high = ci_logit[, 2], odds_ratio = exp(co[, "Estimate"]), OR_CI_low = exp(ci_logit[, 1]), OR_CI_high = exp(ci_logit[, 2]))
write.csv(topo_out, file.path(out_dir, "mpr_logistic_regression.csv"), row.names = FALSE)
ll_full <- as.numeric(logLik(topo_fit))
ll_null <- as.numeric(logLik(glm(PresenceAb ~ 1, family = binomial(), data = mpr)))
mcfadden_r2 <- 1 - ll_full / ll_null

resid_pearson <- residuals(topo_fit, type = "pearson")
resid_neighbors <- knn_indices(mpr$Long, mpr$Lat, k = 8)
resid_i <- moran_i_knn(resid_pearson, resid_neighbors)
set.seed(20261101)
perm <- replicate(9999, moran_i_knn(sample(resid_pearson), resid_neighbors))
resid_p <- (1 + sum(abs(perm) >= abs(resid_i))) / (length(perm) + 1)

# Exact-coordinate overlap accounting.
book_coord <- paste(b09$Long, b09$Lat, sep = "_")
mpr_coord <- paste(mpr$Long, mpr$Lat, sep = "_")
overlap <- data.frame(Book2009_n = nrow(b09), MPR_n = nrow(mpr), exact_coordinate_overlap = length(intersect(book_coord, mpr_coord)),
                      overlap_Book_pct = length(intersect(book_coord, mpr_coord)) / nrow(b09), overlap_MPR_pct = length(intersect(book_coord, mpr_coord)) / nrow(mpr))
write.csv(overlap, file.path(out_dir, "mpr_book_coordinate_overlap.csv"), row.names = FALSE)

# Human-readable analysis summary for manuscript drafting and reviewer response.
summary_lines <- c(
  "REPRODUCIBLE ANALYSIS SUMMARY",
  "",
  "Design: A new, documented matched spatial validation was run for the 2005 and 2009 Book datasets. It is not represented as a reproduction of the undocumented historical spatial analysis.",
  sprintf("Coordinate geometry: %d-by-%d coordinate-unit grid blocks with a %d-unit test-to-training buffer and grid-origin offsets of (%d, %d) coordinate units. The CSVs do not contain a CRS field, so metric units must be confirmed by the authors before publication.", block_size, block_size, buffer_size, outer_shift_x, outer_shift_y),
  "Models: RF used 500 trees; mtry was selected within the outer training data from {2,4,7} using 3-fold inner spatial CV. SVM used e1071 C-classification with a radial kernel, scale=TRUE, no class weights, and cost {0.25,1,4,16} x gamma {1/14,1/7,2/7} selected by inner spatial CV. Epsilon was not used.",
  "No response-informed feature selection or threshold tuning was performed; both models used the same seven archived binary-comparison predictors.",
  "",
  "ARCHIVED RANDOM SPLIT (separate from new spatial CV):",
  paste(capture.output(print(archived_counts, row.names = FALSE)), collapse = "\n"),
  "",
  "NEW MATCHED SPATIAL-CV PERFORMANCE:",
  paste(capture.output(print(perf, row.names = FALSE, digits = 5)), collapse = "\n"),
  "",
  "BLOCK-BOOTSTRAP INTERVALS:",
  paste(capture.output(print(boot, row.names = FALSE, digits = 5)), collapse = "\n"),
  "",
  "FOLD ACCOUNTING:",
  paste(capture.output(print(folds, row.names = FALSE)), collapse = "\n"),
  "",
  "MPR TOPOGRAPHIC ANALYSIS:",
  sprintf("Pearson chi-square = %.5f, df = %d, asymptotic p = %.5g; Monte Carlo p = %.5g (99,999 simulations).", chisq_asymptotic$statistic, chisq_asymptotic$parameter, chisq_asymptotic$p.value, chisq_montecarlo$p.value),
  sprintf("McFadden R2 = %.5f. Pearson-residual Moran's I (8-nearest-neighbor, row-standardized) = %.5f; two-sided permutation p = %.5g (9,999 permutations).", mcfadden_r2, resid_i, resid_p),
  paste(capture.output(print(aspect_summary, row.names = FALSE, digits = 5)), collapse = "\n"),
  paste(capture.output(print(topo_out, row.names = FALSE, digits = 5)), collapse = "\n"),
  "",
  "MPR / BOOK OVERLAP:",
  paste(capture.output(print(overlap, row.names = FALSE, digits = 5)), collapse = "\n")
)
writeLines(summary_lines, file.path(out_dir, "revised_analysis_summary.txt"))

# 5. TRAINING-ONLY COLLINEARITY AND REFERENCE-LABEL MORAN DIAGNOSTICS ------------
# Reviewer-3 diagnostics only: no feature exclusion, RF/SVM refit, or lag predictor.
# Inputs: Book data.frames (or CSV paths), and seed-123 split labels 1=train, 2=holdout.
# Correlations/VIF use training OIDs only; label Moran uses all reference records.
diagnose_knn <- function(x, y, k = 8L) {
  n <- length(x)
  stopifnot(length(y) == n, all(is.finite(x)), all(is.finite(y)),
            length(k) == 1L, is.finite(k), k == as.integer(k), k >= 1L, k < n)
  d2 <- outer(x, x, "-")^2 + outer(y, y, "-")^2
  diag(d2) <- Inf                 # Exclude self even with duplicate coordinates.
  indices <- vapply(seq_len(n), function(i) order(d2[i, ], seq_len(n))[seq_len(k)], integer(k))
  t(matrix(indices, nrow = k, ncol = n))
}

diagnose_moran <- function(z, neighbors) {
  stopifnot(is.matrix(neighbors), nrow(neighbors) == length(z), all(is.finite(z)),
            ncol(neighbors) > 0L, all(neighbors >= 1L & neighbors <= length(z)))
  zz <- z - mean(z)
  denom <- sum(zz^2)
  if (denom <= 0) return(NA_real_)
  lag <- rowMeans(matrix(zz[neighbors], nrow = nrow(neighbors), ncol = ncol(neighbors)))
  sum(zz * lag) / denom            # Row-standardized weights: S0=N, hence N/S0=1.
}

run_reference_diagnostics <- function(b05, b09, split_id, out_dir) {
  as_book <- function(b) if (is.character(b) && length(b) == 1L)
    read.csv(b, check.names = FALSE, stringsAsFactors = FALSE) else b
  b05 <- as_book(b05); b09 <- as_book(b09)
  active <- c("band4", "NBR", "B54", "TCG", "DEM", "AnomB54", "DiffB54")
  candidate <- c("band4", "band5", "NBR", "B54", "TCG", "MeanB54", "DEM", "AnomB54", "DiffB54")
  sets <- list(active7 = active, diagnostic9 = candidate)
  required <- c("OID_", "year", "Long", "Lat", "PresenceAb", "Classcode", "Classname", "Classvalue", candidate)
  for (b in list(b05, b09)) {
    stopifnot(is.data.frame(b), all(required %in% names(b)), !anyNA(b$OID_),
              !anyDuplicated(b$OID_), !anyNA(b$PresenceAb), all(b$PresenceAb %in% 0:1))
  }
  stopifnot(all(b05$year == 2005), all(b09$year == 2009), nrow(b05) == nrow(b09),
            setequal(b05$OID_, b09$OID_), length(split_id) == nrow(b05),
            !anyNA(split_id), all(split_id %in% 1:2), any(split_id == 1L))
  for (nm in c("Long", "Lat", "PresenceAb"))
    stopifnot(isTRUE(all.equal(b05[[nm]], b09[[nm]][match(b05$OID_, b09$OID_)], check.attributes = FALSE)))
  train_oids <- b05$OID_[split_id == 1L]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stopifnot(dir.exists(out_dir))
  save_csv <- function(x, name) write.csv(x, file.path(out_dir, name), row.names = FALSE, na = "NA")
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  bind <- function(x) { r <- do.call(rbind, x); rownames(r) <- NULL; r }
  # Preserve caller RNG state; R permutation stream is not NumPy-bit-identical.
  old_kind <- RNGkind(); had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({ do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  pairs <- summaries <- strongest <- vifs <- scaling <- audits <- morans <- accounting <- crosswalk <- membership <- list()
  for (b in list(b05, b09)) {
    year <- unique(b$year); tr <- b[match(train_oids, b$OID_), , drop = FALSE]
    partition <- ifelse(split_id[match(b$OID_, b05$OID_)] == 1L, "training", "holdout")
    membership[[length(membership) + 1L]] <- data.frame(year, OID_ = b$OID_, partition)
    x <- as.data.frame(lapply(tr[candidate], num), check.names = FALSE)
    for (c in candidate) {
      finite <- is.finite(x[[c]]); u <- length(unique(x[[c]][finite]))
      audits[[length(audits) + 1L]] <- data.frame(year, predictor = c, training_n = nrow(tr),
        raw_missing_n = sum(is.na(tr[[c]])), non_numeric_n = sum(!is.na(tr[[c]]) & is.na(x[[c]])),
        nonfinite_nonmissing_n = sum(!is.na(x[[c]]) & !finite), finite_n = sum(finite),
        finite_unique_n = u, constant = u <= 1L)
      x[[c]][!finite] <- NA_real_
    }
    year_pairs <- list()
    for (set_name in names(sets)) for (method in c("pearson", "spearman")) {
      cols <- sets[[set_name]]; mat <- matrix(NA_real_, length(cols), length(cols), dimnames = list(cols, cols))
      for (c in cols) if (sum(!is.na(x[[c]])) >= 2L && length(unique(na.omit(x[[c]]))) > 1L) mat[c, c] <- 1
      set_pairs <- list()
      combinations <- combn(cols, 2L)
      for (j in seq_len(ncol(combinations))) {
        ab <- combinations[, j]; xy <- x[complete.cases(x[ab]), ab, drop = FALSE]; n <- nrow(xy)
        status <- if (n < 2L) "insufficient_pair_observations" else if
          (any(vapply(xy, function(v) length(unique(v)) <= 1L, logical(1)))) "constant_in_pair" else "ok"
        value <- if (status == "ok") cor(xy[[1L]], xy[[2L]], method = method) else NA_real_
        mat[ab[1L], ab[2L]] <- mat[ab[2L], ab[1L]] <- value
        set_pairs[[j]] <- data.frame(year, predictor_set = set_name, method,
          predictor_1 = ab[1L], predictor_2 = ab[2L], n_pair = n, correlation = value,
          abs_correlation = abs(value), status)
      }
      save_csv(data.frame(predictor = cols, mat, check.names = FALSE), sprintf("correlation_%s_%s_%s.csv", year, set_name, method))
      sp <- bind(set_pairs); valid <- sp[sp$status == "ok", , drop = FALSE]
      valid <- valid[order(-valid$abs_correlation, valid$predictor_1, valid$predictor_2, method = "radix"), , drop = FALSE]
      summaries[[length(summaries) + 1L]] <- data.frame(year, predictor_set = set_name, method, training_n = nrow(tr),
        total_unique_pairs = nrow(sp), valid_unique_pairs = nrow(valid), undefined_pairs = nrow(sp) - nrow(valid),
        abs_r_gt_0_6 = sum(valid$abs_correlation > .6), abs_r_gt_0_7 = sum(valid$abs_correlation > .7),
        abs_r_gt_0_8 = sum(valid$abs_correlation > .8), strongest_pair = if (nrow(valid))
          paste(valid$predictor_1[1L], valid$predictor_2[1L], sep = " / ") else "",
        strongest_correlation = if (nrow(valid)) valid$correlation[1L] else NA_real_)
      pairs[[length(pairs) + 1L]] <- sp; year_pairs[[length(year_pairs) + 1L]] <- sp
    }
    yp <- bind(year_pairs); a <- yp[yp$predictor_set == "active7" & yp$method == "pearson", ]
    s <- yp[yp$predictor_set == "active7" & yp$method == "spearman", ]
    top <- data.frame(year, rank = seq_len(nrow(a)), a[c("predictor_1", "predictor_2", "n_pair")],
                      pearson = a$correlation, spearman = s$correlation)
    top$max_abs_correlation <- pmax(abs(top$pearson), abs(top$spearman), na.rm = TRUE)
    top$max_abs_correlation[is.na(top$pearson) & is.na(top$spearman)] <- NA_real_
    top <- head(top[order(-top$max_abs_correlation, top$predictor_1, top$predictor_2, method = "radix"), ], 10L)
    top$rank <- seq_len(nrow(top)); strongest[[length(strongest) + 1L]] <- top
    xv <- as.matrix(x[complete.cases(x[active]), active, drop = FALSE])
    means <- colMeans(xv); sds <- apply(xv, 2L, sd); enough <- nrow(xv) > length(active) && all(is.finite(sds) & sds > 0)
    if (enough) {
      z <- sweep(sweep(xv, 2L, means, "-"), 2L, sds, "/"); full <- cbind(1, z)
      sing <- svd(full, nu = 0L, nv = 0L)$d; design_rank <- sum(sing > max(dim(full)) * max(sing) * .Machine$double.eps)
    }
    for (c in active) {
      scaling[[length(scaling) + 1L]] <- data.frame(year, predictor = c, n_complete_case = nrow(xv),
                                                 training_mean = means[c], training_sample_sd = sds[c])
      fraction <- if (enough) sum(lm.fit(cbind(1, z[, setdiff(active, c), drop = FALSE]), z[, c])$residuals^2) /
        sum((z[, c] - mean(z[, c]))^2) else NA_real_
      vif <- if (!enough) NA_real_ else if (fraction <= .Machine$double.eps) Inf else 1 / fraction
      vifs[[length(vifs) + 1L]] <- data.frame(year, predictor = c, n_complete_case = nrow(xv),
        R2_auxiliary = 1 - fraction, VIF = vif, design_rank = if (enough) design_rank else NA_integer_,
        status = if (!enough) "not_computed_constant_or_insufficient_complete_cases" else if (is.infinite(vif)) "exact_linear_dependence" else "ok")
    }
    # Source Long/Lat values are projected-coordinate units, NOT longitude/latitude degrees.
    neighbors <- diagnose_knn(num(b$Long), num(b$Lat), 8L); y <- num(b$PresenceAb)
    observed <- diagnose_moran(y, neighbors); set.seed(20261104)
    extreme <- if (is.finite(observed)) sum(vapply(seq_len(9999L), function(i)
      abs(diagnose_moran(sample(y), neighbors)) >= abs(observed), logical(1))) else NA_integer_
    morans[[length(morans) + 1L]] <- data.frame(year, N = nrow(b), k = 8L,
      weights = "directed 8-nearest-neighbor row-standardized", Moran_I = observed,
      two_sided_permutation_p = (extreme + 1) / 10000, permutations = 9999L, seed = 20261104L)
    positive <- sum(y == 1); negative <- sum(y == 0); prevalence <- positive / length(y)
    zcrit <- qnorm(.975); nominal <- zcrit^2 * prevalence * (1 - prevalence) / .05^2
    accounting[[length(accounting) + 1L]] <- data.frame(year, N = nrow(b), positive_n = positive, negative_n = negative,
      training_n = sum(partition == "training"), training_positive_n = sum(y == 1 & partition == "training"),
      training_negative_n = sum(y == 0 & partition == "training"), holdout_n = sum(partition == "holdout"),
      holdout_positive_n = sum(y == 1 & partition == "holdout"), holdout_negative_n = sum(y == 0 & partition == "holdout"),
      prevalence, confidence_level = .95, margin_of_error = .05, z_critical = zcrit,
      cochran_n_unrounded = nominal, cochran_n_ceiling = ceiling(nominal))
    # Record-level crosswalk preserves IDs/raw coding; absent Classvar is never invented.
    crosswalk[[length(crosswalk) + 1L]] <- data.frame(year, OID_ = b$OID_, Classcode = b$Classcode,
      Classname = b$Classname, Classvalue = b$Classvalue,
      Classvar = if ("Classvar" %in% names(b)) b$Classvar else NA, PresenceAb = b$PresenceAb,
      partition, Classvar_present = "Classvar" %in% names(b), check.names = FALSE)
  }
  tidy <- bind(pairs); tidy <- tidy[order(tidy$year, tidy$predictor_set, tidy$method, -tidy$abs_correlation, method = "radix"), ]
  outputs <- list(collinearity_tidy_pairs = tidy, collinearity_threshold_summary = bind(summaries),
    collinearity_strongest_active7_pairs = bind(strongest), collinearity_vif_active7 = bind(vifs),
    collinearity_vif_training_scaling = bind(scaling), collinearity_column_audit = bind(audits),
    reference_label_moran = bind(morans), reference_sample_accounting = bind(accounting),
    reference_code_crosswalk = bind(crosswalk), reference_diagnostic_membership = bind(membership))
  for (name in names(outputs)) save_csv(outputs[[name]], paste0(name, ".csv"))
  invisible(outputs)
}

run_reference_diagnostics(b05, b09, archived_split, out_dir)

# 6. OPTIONAL GEOMETRY SENSITIVITY (TABLE S5) -----------------------------------
# Primary geometry above always runs. Set RUN_SENSITIVITY=TRUE to rerun the three
# additional geometries actually included in the manuscript. No synthetic
# data, repeated random allocations, or new feature-removal experiment is used.
extract_sensitivity <- function(label, perf, boot) {
  do.call(rbind, lapply(c("2005", "2009"), function(y) {
    p <- perf[perf$year == y, ]; b <- boot[boot$year == y & boot$metric == "Delta_MCC_RF_minus_SVM", ]
    data.frame(geometry = label, year = y,
               RF_MCC = p$MCC[p$model == "RF"], SVM_MCC = p$MCC[p$model == "SVM"],
               delta_RF_minus_SVM = b$estimate, delta_CI_low = b$CI_low, delta_CI_high = b$CI_high)
  }))
}
sensitivity <- list(extract_sensitivity("5k block / 1k buffer, origin (0,0)", perf, boot))
if (RUN_SENSITIVITY) {
  geometries <- data.frame(
    label = c("2.5k block / 0.5k buffer", "5k block / 1k buffer, origin (0,2.5k)", "7.5k block / 1.5k buffer"),
    block = c(2500, 5000, 7500), buffer = c(500, 1000, 1500),
    shift_x = c(0, 0, 0), shift_y = c(0, 2500, 0))
  for (g in seq_len(nrow(geometries))) {
    block_size <- geometries$block[g]; buffer_size <- geometries$buffer[g]
    outer_shift_x <- geometries$shift_x[g]; outer_shift_y <- geometries$shift_y[g]
    message("Sensitivity geometry: ", geometries$label[g])
    result05 <- run_spatial_nested(b05, "2005", 1)
    result09 <- run_spatial_nested(b09, "2009", 2)
    oo <- rbind(result05$oof, result09$oof)
    pp <- do.call(rbind, lapply(split(oo, oo$year), function(x) {
      data.frame(year = unique(x$year), model = c("RF", "SVM"),
                 rbind(metric_row(x$observed, x$RF_prediction),
                       metric_row(x$observed, x$SVM_prediction)))
    }))
    bb <- do.call(rbind, lapply(split(oo, oo$year), function(x) {
      cbind(year = unique(x$year), block_bootstrap(x, x$outer_block, B = 9999,
                                                  seed = 450 + as.integer(unique(x$year))))
    }))
    sensitivity[[length(sensitivity) + 1]] <- extract_sensitivity(geometries$label[g], pp, bb)
    folder <- file.path(out_dir, paste0("geometry_", g))
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    write.csv(oo, file.path(folder, "spatial_oof_predictions.csv"), row.names = FALSE)
    write.csv(pp, file.path(folder, "spatial_performance_metrics.csv"), row.names = FALSE)
    write.csv(bb, file.path(folder, "spatial_block_bootstrap.csv"), row.names = FALSE)
    write.csv(rbind(result05$folds, result09$folds), file.path(folder, "spatial_fold_accounting.csv"), row.names = FALSE)
    write.csv(rbind(result05$tuning, result09$tuning), file.path(folder, "nested_tuning_results.csv"), row.names = FALSE)
  }
}
sens <- do.call(rbind, sensitivity)
filename <- if (RUN_SENSITIVITY) "spatial_geometry_sensitivity.csv" else "spatial_geometry_primary_only.csv"
write.csv(sens, file.path(out_dir, filename), row.names = FALSE)
if (!RUN_SENSITIVITY) message("Only the primary geometry was run. Set RUN_SENSITIVITY <- TRUE for full Table S5.")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
write.csv(data.frame(file = basename(files), MD5 = unname(tools::md5sum(files))),
          file.path(out_dir, "input_file_checksums.csv"), row.names = FALSE)
message("Finished. Results are in: ", normalizePath(out_dir))
