# Reproducible analysis prepared in response to geosci-737 review
# Inputs are unmodified working copies of supplied CSV files.
# This script intentionally distinguishes a new, documented analysis from historical results
# whose SVM code, spatial blocks, and prediction records were not retained.

suppressPackageStartupMessages({
  library(randomForest)
  library(e1071)
})

root <- "/home/ubuntu/revision_work"
input_dir <- file.path(root, "input")
out_dir <- Sys.getenv("OUTPUT_DIR", unset = file.path(root, "analysis"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

features <- c("band4", "NBR", "B54", "TCG", "DEM", "AnomB54", "DiffB54")
block_size <- as.numeric(Sys.getenv("BLOCK_SIZE", unset = "5000")) # coordinate units; intended as metres only after CRS confirmation
buffer_size <- as.numeric(Sys.getenv("BUFFER_SIZE", unset = "1000")) # coordinate units; intended as metres only after CRS confirmation
outer_shift_x <- as.numeric(Sys.getenv("OUTER_SHIFT_X", unset = "0"))
outer_shift_y <- as.numeric(Sys.getenv("OUTER_SHIFT_Y", unset = "0"))
outer_k <- 5
inner_k <- 3
rf_mtry_grid <- c(2, 4, 7)
svm_cost_grid <- c(0.25, 1, 4, 16)
svm_gamma_grid <- c(1/14, 1/7, 2/7)

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
    best <- grid[which.max(replace(grid$inner_MCC, is.na(grid$inner_MCC), -Inf)), , drop = FALSE]
  } else {
    grid <- expand.grid(cost = svm_cost_grid, gamma = svm_gamma_grid)
    grid$inner_MCC <- vapply(seq_len(nrow(grid)), function(i) {
      inner_score(train, "SVM", grid[i, , drop = FALSE], year_index, outer_fold)
    }, numeric(1))
    best <- grid[which.max(replace(grid$inner_MCC, is.na(grid$inner_MCC), -Inf)), , drop = FALSE]
  }
  list(grid = grid, best = best)
}

run_spatial_nested <- function(dat, year_label, year_index) {
  dat$block_outer <- make_blocks(dat, block_size, shift_x = outer_shift_x, shift_y = outer_shift_y)
  dat$fold_outer <- assign_folds(dat$block_outer, dat$PresenceAb, outer_k)
  oof <- list(); tuning <- list(); fold_info <- list()
  for (f in seq_len(outer_k)) {
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

# Load binary classification datasets
b05 <- read.csv(file.path(input_dir, "Book2005.csv"), check.names = FALSE)
b09 <- read.csv(file.path(input_dir, "Book2009.csv"), check.names = FALSE)
stopifnot(all(features %in% names(b05)), all(features %in% names(b09)))
stopifnot(identical(b05$OID_, b09$OID_), identical(b05$Long, b09$Long), identical(b05$Lat, b09$Lat), identical(b05$PresenceAb, b09$PresenceAb))

# Document archived random split, which is separate from the new spatial-CV results.
set.seed(123)
archived_split <- sample(2, nrow(b05), replace = TRUE, prob = c(0.7, 0.3))
archived_counts <- data.frame(partition = c("training", "test"),
                              n = c(sum(archived_split == 1), sum(archived_split == 2)),
                              negative = c(sum(b05$PresenceAb[archived_split == 1] == 0), sum(b05$PresenceAb[archived_split == 2] == 0)),
                              positive = c(sum(b05$PresenceAb[archived_split == 1] == 1), sum(b05$PresenceAb[archived_split == 2] == 1)))
write.csv(archived_counts, file.path(out_dir, "archived_seed123_split_counts.csv"), row.names = FALSE)

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

# MPR topographic analysis, reported explicitly as a supplemental overlapping sample.
mpr <- read.csv(file.path(input_dir, "MPR_all2009.csv"), check.names = FALSE)
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

# Residual Moran's I with explicitly specified 8-nearest-neighbor, row-standardized weights.
# The neighbour matrix is calculated once; permutations reuse it rather than materializing dense matrices.
knn_indices <- function(x, y, k = 8) {
  n <- length(x); d2 <- outer(x, x, "-")^2 + outer(y, y, "-")^2; diag(d2) <- Inf
  t(apply(d2, 1, function(z) order(z)[seq_len(k)]))
}
moran_i_knn <- function(z, neighbors) {
  zz <- z - mean(z)
  sum(zz * rowMeans(matrix(zz[neighbors], nrow = nrow(neighbors), ncol = ncol(neighbors)))) / sum(zz^2)
}
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
  "No response-informed feature selection or threshold tuning was performed; both models used the same seven archived-RF predictors.",
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
