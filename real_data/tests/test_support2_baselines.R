#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
root <- dirname(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(quantreg))
source(file.path(root, "utils/utils.R"))
source(file.path(root, "utils/cqr.R"))
source(file.path(root, "real_data/support2.R"))
source(file.path(root, "real_data/baselines.R"))
check <- function(label, value) {
  if (!isTRUE(value)) stop("FAILED: ", label)
  cat("PASS: ", label, "\n", sep = "")
}

# Independent comparison against the historical calibration formulas on a
# controlled event learner and the same explicit candidate grid.
predict.support2_fixture <- function(object, newdata, type = "quantile", p = NULL, ...) {
  if (object$kind == "censor") return(log(12) + .3 * newdata$x)
  if (object$kind == "forest") return(list(predictions = matrix(
    exp(log(12) + .3 * newdata$x + .5 * qnorm(.95)), ncol = 1)))
  outer(40 + 5 * newdata$x, p, "*")
}
local({
  mdl <- structure(list(kind = "event"), class = "support2_fixture")
  sg <- structure(list(fit = structure(list(kind = "censor", scale = .5),
    class = "support2_fixture"), xnames = "x"), class = "sc_aft_lognormal")
  cal <- data.frame(x = seq(-.5, .5, length.out = 40),
    censored_T = seq(2, 35, length.out = 40), event = rep(c(1L, 0L, 1L, 1L), 10))
  cal$C <- cal$censored_T + ifelse(cal$event == 1, 20, 0)
  test <- data.frame(x = c(-.4, 0, .4))
  grid <- seq(.02, .98, .02)
  env <- new.env(parent = .GlobalEnv)
  sys.source(file.path(root, "utils/est_alpha.R"), env)
  env$v_pts_qt <- env$v_pts_qct <- function(...) 1 - grid
  env$sc_prob <- function(...) pmax(sc_prob(...), .05)
  env$lv_qt <- function(mdl, newdata, v, ...) as.numeric(predict(mdl, newdata, p = 1 - v))
  env$lv_qct <- function(mdl, newdata, v, alpha, qc) pmin(
    as.numeric(predict(mdl, newdata, p = 1 - v)), qc)
  old <- env$alpha_qt(mdl, test, cal, cal, "x", .1, 3, sg, .95)$lower_bnd_l
  new <- support2_dft_adaptive(mdl, sg, cal, test, cal$C, "x", .1, .05, grid, 180)$bound
  check("adaptive DFT-T matches historical risk-ratio calibration", max(abs(old - new)) < 1e-10)
  forest <- structure(list(kind = "forest"), class = "support2_fixture")
  old_ct <- env$alpha_qct(mdl, forest, test, cal, cal, "x", .1, 3, sg, .95)$lower_bnd_l
  new_ct <- support2_dft_adaptive(mdl, sg, cal, test, cal$C, "x", .1, .05, grid, 180, TRUE)$bound
  check("adaptive DFT-CT matches historical capped calibration", max(abs(old_ct - new_ct)) < 1e-10)
  none <- cal; none$censored_T <- .001; none$C <- .002; none$event <- 1L
  fallback <- support2_dft_adaptive(mdl, sg, none, test, none$C, "x", .1, .05, grid, 180)
  check("no eligible adaptive candidate is a recorded zero fallback", all(fallback$bound == 0) && nzchar(fallback$fallback))
  fixed <- support2_dft_fixed(mdl, sg, cal, test, cal$C, "x", .1, .05, 30, 180)
  keep <- cal$C >= 30
  scores <- pmin(30, as.numeric(predict(mdl, cal[keep, ], p = .1))) - pmin(30, cal$censored_T[keep])
  wc <- 1 / pmax(sc_prob(sg, cal[keep, ], "x", 30), .05)
  wt <- 1 / pmax(sc_prob(sg, test, "x", 30), .05)
  correction <- vapply(wt, function(w) get_calibration(scores, wc, w, .1), numeric(1))
  expected <- pmax(0, pmin(30, as.numeric(predict(mdl, test, p = .1))) - correction)
  check("fixed DFT preserves weighted order statistic and test-point mass", max(abs(fixed$bound - expected)) < 1e-10)
  empty <- support2_dft_fixed(mdl, sg, cal, test, cal$C, "x", .1, .05, 170, 180)
  check("empty fixed-DFT selection returns a documented conservative bound", all(empty$bound == 0) && empty$n_selected == 0 && nzchar(empty$fallback))
})

raw <- load_support2(cache_dir = file.path(root, "real_data/data"))
p <- prepare_support2(raw, seed = 13, scenario = "native", horizon = 730)
mdl <- support2_checked_survreg(p$data_fit, p$xnames)
sg <- support2_checked_survreg(p$data_fit, p$xnames, TRUE)
set.seed(55); old_rng <- .Random.seed
C <- support2_impute_c(p$data_calib, sg, 1234)
check("imputation preserves the caller random stream", identical(.Random.seed, old_rng))
check("conditional imputation is reproducible", identical(C, support2_impute_c(p$data_calib, sg, 1234)))
idx <- p$data_calib$event == 0
check("observed censoring times stay unchanged", identical(C[idx], p$data_calib$censored_T[idx]))
check("imputed C strictly exceeds observed death times", all(C[!idx] > p$data_calib$censored_T[!idx]))
mu <- as.numeric(predict(sg$fit, p$data_calib[!idx, ], type = "lp"))
logGc <- pnorm((log(C[!idx]) - mu) / sg$fit$scale, lower.tail = FALSE, log.p = TRUE)
logGy <- pnorm((log(p$data_calib$censored_T[!idx]) - mu) / sg$fit$scale, lower.tail = FALSE, log.p = TRUE)
check("imputed draws follow the truncated conditional survival distribution", abs(mean(exp(logGc - logGy)) - .5) < .025)
before <- p
b <- support2_baselines(p, mdl, sg, .1, .05, seq(.05, .95, .05),
  fixed_cutoff = "180", imputation_seed = 444)
required <- as.vector(outer(c("Adaptive DFT-T", "Adaptive DFT-CT", "DFT-fixed", "CQR"),
  c("joint", "subgroup"), function(m, s) paste0(m, " (", s, ")")))
check("all eight requested method names are present", setequal(names(b$bounds), required))
check("baseline preparation leaves original records unchanged", identical(p, before) && !"C" %in% names(p$data_fit))
check("imputation excludes test patients", !any(b$imputation_manifest$patient_id %in% p$evaluation$patient_id))
altered <- p
altered$data_test$censored_T <- rev(altered$data_test$censored_T) + 1000
altered$data_test$event <- 1L - altered$data_test$event
altered$evaluation$original_Y <- 1
altered$evaluation$original_event <- 0L
b2 <- support2_baselines(altered, mdl, sg, .1, .05, seq(.05, .95, .05),
  fixed_cutoff = "180", imputation_seed = 444)
check("changing test outcomes cannot change any baseline prediction", identical(b$bounds, b2$bounds))
bad <- support2_capture(support2_impute_c(p$data_calib, NULL, 1))
check("unsupported censoring models fail explicitly", inherits(bad$value, "error"))
p3 <- prepare_support2(raw, seed = 3, scenario = "native", horizon = 730)
older <- p3$data_fit[p3$fit_groups == "female__age_ge65", , drop = FALSE]
sx3 <- support2_subgroup_columns(older, p3$xnames)
unpenalized <- support2_capture(support2_checked_survreg(older, sx3, TRUE, ridge_fallback = FALSE))
check("sparse-cell undefined censoring coefficient is detected", inherits(unpenalized$value, "error"))
stabilized <- support2_checked_survreg(older, sx3, TRUE, ridge_fallback = TRUE)
check("training-only ridge fallback resolves the sparse-cell fit", all(is.finite(coef(stabilized$fit))))
check("ridge fallback is explicit in fitted-model metadata", grepl("ridge fallback", support2_model_note(stabilized)))
check("stabilized model predicts valid censoring survival probabilities", all(is.finite(
  sc_prob(stabilized, older, sx3, 730))))
cat("All baseline tests passed.\n")
