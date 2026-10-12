# Meaningful regression and identification checks, using base R assertions.
# Usage: Rscript real_data/tests/test_support2.R [optional original CSV path]
file_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
root <- dirname(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))))
suppressPackageStartupMessages(library(survival))
source(file.path(root, "utils/utils.R"))
source(file.path(root, "utils/cams.R"))
source(file.path(root, "real_data/support2.R"))
check <- function(label, value) {
  if (!isTRUE(value)) stop("FAILED: ", label)
  cat("PASS: ", label, "\n", sep = "")
}

# The compensator must be integrated even before the first censoring jump.
# For exponential S=exp(-r*t), G=exp(-c*t), the integrals are analytical.
local({
  env <- new.env(parent = environment())
  for (name in c("make_aipcw_dense_cache", "aipcw_dense_pseudo_loss")) {
    env[[name]] <- get(name)
    environment(env[[name]]) <- env
  }
  event_rate <- 0.02; censor_rate <- 0.01
  env$sc_prob <- function(mdl0, data, xnames, t) exp(-censor_rate * t)
  env$augmentation_survival_prob <- function(mdl, newdata, t) exp(-event_rate * t)
  dat <- data.frame(x = 0, censored_T = c(50, 4, 3), event = c(0L, 1L, 0L))
  cache <- env$make_aipcw_dense_cache(NULL, dat, "x",
    structure(list(), class = "sc_aft_lognormal"), 1e-8, 20, 2048)
  L <- c(2, 10, 10)
  actual <- env$aipcw_dense_pseudo_loss(L, cache)
  S_L <- exp(-event_rate * L); upper <- pmin(dat$censored_T, L)
  A <- exp(censor_rate * upper) - 1
  B <- censor_rate / (censor_rate + event_rate) *
    (exp((censor_rate + event_rate) * upper) - 1)
  eta <- ifelse(dat$censored_T < L,
    1 - exp(-event_rate * (L - dat$censored_T)), 0)
  expected <- dat$event * (dat$censored_T < L) * exp(censor_rate * dat$censored_T) +
    (1 - dat$event) * eta * exp(censor_rate * dat$censored_T) - A + S_L * B
  check("continuous AIPCW matches analytic exponential integrals",
    max(abs(actual - expected)) < 1e-7)
  check("zero bound has zero pseudo-loss", all(env$aipcw_dense_pseudo_loss(rep(0, 3), cache) == 0))
  check("augmentation is nonzero below first censoring time", abs(actual[1]) > 1e-5)
})

# Test identification and preprocessing using the complete original cohort.
args <- commandArgs(trailingOnly = TRUE)
csv <- if (length(args)) args[1] else NULL
raw <- load_support2(csv, file.path(root, "real_data/data"))
check("original cohort size and outcomes", nrow(raw) == 9105 && sum(raw$death) == 6201)
check("180-day restricted outcome known for every original patient",
  all(raw$d.time[raw$death == 0] >= 180))
p <- prepare_support2(raw, seed = 17, scenario = "added_group")
check("same split seed is reproducible",
  identical(p$split_manifest$split, prepare_support2(raw, seed = 17)$split_manifest$split))
check("split IDs are unique and exhaust cohort", !anyDuplicated(p$split_manifest$patient_id) &&
  sum(vapply(p[c("data_fit", "data_calib", "data_test")], nrow, integer(1))) == 9105)
check("fitting data has no original outcomes or censoring time",
  !any(c("original_Y", "original_event", "restricted_outcome", "T", "C") %in% names(p$data_fit)))
check("added censoring only shortens observed time",
  all(p$split_manifest$observed_Y <= p$split_manifest$original_Y))
check("added event indicators are coherent",
  all(p$split_manifest$observed_event == 0 |
    (p$split_manifest$original_event == 1 &
     p$split_manifest$observed_Y == p$split_manifest$original_Y)))
bad <- try(prepare_support2(raw, scenario = "added_group", horizon = 730), silent = TRUE)
check("unknown restricted outcomes rejected for artificial censoring", inherits(bad, "try-error"))
bad <- try(prepare_support2(raw, feature_set = "day3"), silent = TRUE)
check("day-3 predictors require a landmark", inherits(bad, "try-error"))
edge <- raw[1:2, ]; edge$age <- c(64.999, 65); edge$sex <- "female"
check("age-65 group boundary is correct", identical(support2_groups(edge),
  c("female__age_lt65", "female__age_ge65")))
training <- data.frame(x = c(1, NA_real_, 3))
encoder <- support2_fit_encoder(training, "x")
encoded <- support2_encode(data.frame(x = c(NA_real_, 1e9)), encoder)
check("imputation learned on fitting data only", encoder$specs[[1]]$fill == 2 && encoded$x[1] == 0)

# Hidden natural-censoring outcomes must not be called exact coverage.
local({
  env <- new.env(parent = environment())
  env$support2_evaluate <- support2_evaluate
  environment(env$support2_evaluate) <- env
  env$make_aipcw_cache <- function(..., data_calib) list(G_y = rep(1, nrow(data_calib)))
  env$aipcw_risk_for_bound <- function(lower_bound, ...) rep(0, length(lower_bound))
  env$sc_prob <- function(mdl0, data, ...) rep(1, nrow(data))
  prep <- list(data_test = data.frame(censored_T = c(5, 15), event = c(0, 1)),
    evaluation = data.frame(patient_id = c("a", "b"), original_Y = c(5, 15),
      original_event = c(0, 1), group = c("g", "g")),
    horizon = 20, xnames = "x")
  result <- list(bounds = data.frame(method = c(10, 10)), event_model = NULL,
    censoring_model = NULL)
  m <- env$support2_evaluate(prep, result, 0.05, 128)$metrics
  check("unknown natural outcomes yield NA exact coverage", all(is.na(m$coverage_exact)))
  check("observable coverage brackets retain unknown observations",
    all(m$coverage_observed_lower == 0.5 & m$coverage_observed_upper == 1 &
        m$unknown_indicator_fraction == 0.5))
  prep$evaluation$group <- c("success", "failed")
  result$bounds$method[2] <- NA_real_
  partial <- env$support2_evaluate(prep, result, 0.05, 128)$metrics
  check("one subgroup failure does not erase other group's scores",
    is.finite(partial$coverage_ipcw[partial$audit_group == "success"]) &&
    is.na(partial$coverage_ipcw[partial$audit_group == "failed"]) &&
    is.na(partial$coverage_ipcw[partial$audit_group == "all"]))
})

# Explicit groups should equal the historical X1/sign-X2 selection under legacy integration.
# A saved copy of the unmodified cams.R can be supplied for direct regression.
set.seed(72)
gen <- function(n) {
  x1 <- rep(0:1, length.out = n)
  x2 <- rnorm(n)
  T <- rweibull(n, shape = 2, scale = exp(3 + 0.2 * x1 + 0.1 * x2))
  C <- rlnorm(n, 4, 0.7)
  data.frame(X1 = x1, X2 = x2, censored_T = pmin(T, C), event = as.integer(T <= C))
}
fit <- gen(220); cal <- gen(180); test <- gen(90)
g <- fit_sc_model("aft_lognormal", fit, c("X1", "X2"))
common <- list(x = test[, c("X1", "X2")], p = 2, len_x = nrow(test),
  xnames = c("X1", "X2"), data_fit = fit, data_calib = cal, mdl0 = g,
  alpha = 0.1, augmentation_method = "same", use_intersectional_R = TRUE)
a <- do.call(cams, common)
b <- do.call(cams, c(common, list(calib_groups = make_cams_R_label(cal, TRUE),
                                 test_groups = make_cams_R_label(test, TRUE))))
check("explicit labels preserve legacy four-group predictions", isTRUE(all.equal(a$output, b$output)))
original_path <- Sys.getenv("CAMS_ORIGINAL_CAMS_R")
if (nzchar(original_path)) {
  old <- new.env(parent = .GlobalEnv)
  sys.source(original_path, envir = old)
  original <- do.call(old$cams, common)
  check("existing simulation defaults match unmodified cams.R",
    isTRUE(all.equal(a$output, original$output)) &&
    isTRUE(all.equal(a$diagnostics, original$diagnostics)))
}
cat("All SUPPORT2 tests passed.\n")
