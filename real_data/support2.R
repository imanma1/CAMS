# SUPPORT2 adapter: reuses utils/cams.R, utils/utils.R and utils/cqr.R.
# Model-facing data contains covariates, censored_T and event, never latent T/C.

support2_source_url <- "https://hbiostat.org/data/repo/support2csv.zip"

load_support2 <- function(csv_path = NULL, cache_dir = "real_data/data") {
  if (is.null(csv_path)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    csv_path <- file.path(cache_dir, "support2.csv")
    if (!file.exists(csv_path)) {
      zip_path <- file.path(cache_dir, "support2csv.zip")
      download.file(support2_source_url, zip_path, mode = "wb", quiet = FALSE)
      unzip(zip_path, files = "support2.csv", exdir = cache_dir)
    }
  }
  if (!file.exists(csv_path)) stop("SUPPORT2 CSV does not exist: ", csv_path)
  # Original CSV has an unnamed patient-ID field before its 47 named columns.
  raw <- read.csv(csv_path, row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE, na.strings = c("", "NA"))
  required <- c("age", "sex", "race", "dzgroup", "num.co", "diabetes",
                "dementia", "d.time", "death")
  if (!all(required %in% names(raw))) stop("Unexpected SUPPORT2 CSV schema.")
  raw$patient_id <- rownames(raw)
  for (name in names(raw)) {
    if (is.character(raw[[name]])) raw[[name]] <- trimws(raw[[name]])
  }
  if (anyNA(raw[, c("age", "sex", "d.time", "death")]) ||
      any(!raw$sex %in% c("female", "male")) ||
      any(!raw$death %in% c(0, 1)) || any(raw$d.time <= 0) ||
      anyDuplicated(raw$patient_id)) stop("Invalid outcome, group attribute or ID.")
  attr(raw, "source_path") <- normalizePath(csv_path)
  attr(raw, "source_md5") <- unname(tools::md5sum(csv_path))
  raw
}

support2_groups <- function(raw, age_cutoff = 65) {
  paste(raw$sex, ifelse(raw$age < age_cutoff, paste0("age_lt", age_cutoff),
                        paste0("age_ge", age_cutoff)), sep = "__")
}

support2_split <- function(groups, seed, proportions = c(0.5, 0.25, 0.25)) {
  if (length(proportions) != 3L || any(proportions <= 0) ||
      abs(sum(proportions) - 1) > 1e-10) stop("Invalid split proportions.")
  set.seed(seed)
  labels <- rep(NA_character_, length(groups))
  for (g in sort(unique(groups))) {
    indices <- sample(which(groups == g))
    n_fit <- floor(length(indices) * proportions[1])
    n_cal <- floor(length(indices) * proportions[2])
    if (min(n_fit, n_cal, length(indices) - n_fit - n_cal) < 2L) {
      stop("Too few observations in a split for group ", g)
    }
    labels[indices] <- c(rep("fit", n_fit), rep("calibration", n_cal),
                         rep("test", length(indices) - n_fit - n_cal))
  }
  labels
}

# Learn every fill value, category vocabulary and numerical scale on fit only.
support2_fit_encoder <- function(raw_fit, predictors) {
  specs <- lapply(predictors, function(name) {
    x <- raw_fit[[name]]
    if (is.numeric(x)) {
      finite <- is.finite(x)
      if (!any(finite)) stop("No fitting values for predictor ", name)
      fill <- median(x[finite])
      x[!finite] <- fill
      scale <- sd(x)
      if (!is.finite(scale) || scale == 0) scale <- 1
      list(name = name, type = "numeric", fill = fill, center = mean(x),
           scale = scale, add_missing = any(!finite))
    } else {
      x <- as.character(x)
      missing <- is.na(x) | !nzchar(x)
      available <- x[!missing]
      if (!length(available)) stop("No fitting values for predictor ", name)
      reference <- names(sort(table(available), decreasing = TRUE))[1]
      x[missing] <- ".missing"
      levels <- c(reference, sort(setdiff(unique(x), reference)))
      list(name = name, type = "categorical", levels = levels,
           reference = reference)
    }
  })
  encoder <- list(specs = specs, columns = NULL)
  X <- support2_encode(raw_fit, encoder)
  design <- cbind(.intercept = 1, as.matrix(X))
  q <- qr(design)
  keep <- q$pivot[seq_len(q$rank)]
  encoder$columns <- setdiff(colnames(design)[keep], ".intercept")
  encoder$dropped_columns <- setdiff(names(X), encoder$columns)
  if (!length(encoder$columns)) stop("No estimable predictors.")
  encoder
}

support2_encode <- function(raw, encoder) {
  columns <- list()
  for (spec in encoder$specs) {
    x <- raw[[spec$name]]
    if (spec$type == "numeric") {
      missing <- !is.finite(x)
      x[missing] <- spec$fill
      columns[[spec$name]] <- (x - spec$center) / spec$scale
      if (spec$add_missing) columns[[paste0(spec$name, "__missing")]] <-
        as.numeric(missing)
    } else {
      x <- as.character(x)
      x[is.na(x) | !nzchar(x)] <- ".missing"
      # Categories absent from fitting data use its learned reference category.
      x[!x %in% spec$levels] <- spec$reference
      for (level in spec$levels[-1L]) {
        columns[[paste0(spec$name, "__", level)]] <- as.numeric(x == level)
      }
    }
  }
  X <- as.data.frame(columns, check.names = FALSE)
  names(X) <- make.names(names(X), unique = TRUE)
  if (!is.null(encoder$columns)) X <- X[, encoder$columns, drop = FALSE]
  X
}

prepare_support2 <- function(raw, seed = 1L, scenario = "native", horizon = 180,
                             feature_set = "entry", landmark_days = 0,
                             age_cutoff = 65) {
  scenario <- match.arg(scenario, c("native", "added_common", "added_group"))
  feature_set <- match.arg(feature_set, c("entry", "day3"))
  if (!is.finite(horizon) || horizon <= 0 || !is.finite(landmark_days) ||
      landmark_days < 0) stop("Invalid horizon or landmark.")
  if (feature_set == "day3" && landmark_days < 3) {
    stop("Day-3 physiology requires an explicit landmark of at least three days.")
  }
  cohort <- raw[raw$d.time > landmark_days, , drop = FALSE]
  Y_original <- cohort$d.time - landmark_days
  D_original <- as.integer(cohort$death)
  groups <- support2_groups(cohort, age_cutoff)
  split <- support2_split(groups, seed)
  predictors <- c("age", "sex", "race", "dzgroup", "num.co", "diabetes", "dementia")
  if (feature_set == "day3") predictors <- c(predictors, "scoma", "meanbp", "wblc",
    "hrt", "resp", "temp", "pafi", "alb", "bili", "crea", "sod", "ph")
  if (!all(predictors %in% names(cohort))) stop("Missing required predictors.")
  encoder <- support2_fit_encoder(cohort[split == "fit", , drop = FALSE], predictors)
  X <- support2_encode(cohort, encoder)
  Y <- Y_original
  D <- D_original
  A <- rep(NA_real_, length(Y))
  if (scenario != "native") {
    # Check the original cohort, not just those assigned to the test set.
    if (any(D_original == 0 & Y_original < horizon)) {
      stop("Exact added-censoring experiment requires original follow-up through horizon.")
    }
    rho <- if (scenario == "added_common") rep(0.8, length(Y)) else {
      ifelse(cohort$sex == "female", ifelse(cohort$age < age_cutoff, 0.8, 0.65),
              ifelse(cohort$age < age_cutoff, 0.65, 0.5))
    }
    set.seed(as.integer(seed) + if (scenario == "added_common") 100001L else 200001L)
    A <- rexp(length(Y), rate = -log(rho) / horizon)
    Y0 <- pmin(Y_original, horizon)
    D0 <- as.integer(D_original == 1 & Y_original <= horizon)
    Y <- pmin(Y0, A)
    D <- as.integer(D0 == 1 & Y0 <= A)
  }
  data <- cbind(X, censored_T = Y, event = D)
  parts <- lapply(c("fit", "calibration", "test"), function(p) {
    data[split == p, , drop = FALSE]
  })
  names(parts) <- c("data_fit", "data_calib", "data_test")
  # This table is used only after predictions are fixed.
  evaluation <- data.frame(patient_id = cohort$patient_id[split == "test"],
    original_Y = Y_original[split == "test"],
    original_event = D_original[split == "test"],
    restricted_outcome = pmin(Y_original[split == "test"], horizon),
    group = groups[split == "test"], stringsAsFactors = FALSE)
  manifest <- data.frame(patient_id = cohort$patient_id, split = split, group = groups,
    original_Y = Y_original, original_event = D_original, observed_Y = Y,
    observed_event = D, added_censoring_time = A, stringsAsFactors = FALSE)
  c(parts, list(xnames = names(X), calib_groups = groups[split == "calibration"],
    test_groups = groups[split == "test"], fit_groups = groups[split == "fit"],
    evaluation = evaluation, split_manifest = manifest, encoder = encoder,
    predictors = predictors, n_excluded_landmark = nrow(raw) - nrow(cohort),
    scenario = scenario, horizon = horizon))
}

support2_model_formula <- function(xnames) {
  as.formula(paste("Surv(censored_T, event) ~", paste(xnames, collapse = " + ")))
}

support2_subgroup_columns <- function(data, xnames) {
  design <- cbind(.intercept = 1, as.matrix(data[, xnames, drop = FALSE]))
  q <- qr(design)
  setdiff(colnames(design)[q$pivot[seq_len(q$rank)]], ".intercept")
}

support2_methods <- function(prepared, alpha, probability_floor, candidate_grid,
                             integration_steps, include_subgroups = TRUE,
                             include_cqr = TRUE, include_dft = TRUE,
                             dft_fixed_cutoff = "auto", imputation_seed = 300001L) {
  fit <- prepared$data_fit; cal <- prepared$data_calib; test <- prepared$data_test
  xnames <- prepared$xnames; cap <- prepared$horizon
  start <- proc.time()[3]
  mdl0 <- support2_checked_survreg(fit, xnames, censoring = TRUE)
  mdl <- support2_checked_survreg(fit, xnames)
  if (any(!is.finite(coef(mdl))) || any(!is.finite(coef(mdl0$fit)))) {
    stop("Nonfinite pooled model coefficients.")
  }
  result <- cams(x = test[, xnames, drop = FALSE], p = length(xnames),
    len_x = nrow(test), xnames = xnames, data_fit = fit, data_calib = cal,
    mdl0 = mdl0, alpha = alpha, augmentation_method = "same",
    event_model = mdl, calib_groups = prepared$calib_groups,
    test_groups = prepared$test_groups, bound_cap = cap,
    probability_floor = probability_floor, candidate_grid = candidate_grid,
    integration_method = "dense", integration_steps = integration_steps)
  marginal <- est_alpha_ipcw(mdl, test[, xnames, drop = FALSE], cal, xnames,
    alpha, nrow(test), mdl0, probability_floor, augmentation_mdl = mdl,
    bound_cap = cap, candidate_grid = candidate_grid, integration_method = "dense",
    integration_steps = integration_steps)
  bounds <- data.frame(
    "Weibull-uncalibrated" = pmin(cap, pmax(0, as.numeric(lv_cams(mdl, test, alpha)))),
    "Marginal-IPCW" = marginal$bounds$CAMS,
    "Marginal-AIPCW" = marginal$bounds[["AIPCW-CAMS"]],
    "CAMS-IPCW" = result$output$CAMS,
    "CAMS-AIPCW" = result$output[["AIPCW-CAMS"]], check.names = FALSE)
  diagnostics <- result$diagnostics[result$diagnostics$method %in%
    c("CAMS", "AIPCW-CAMS"), , drop = FALSE]
  diagnostics$method <- ifelse(diagnostics$method == "CAMS", "CAMS-IPCW", "CAMS-AIPCW")
  md <- marginal$diagnostics[marginal$diagnostics$method %in%
    c("CAMS", "AIPCW-CAMS"), , drop = FALSE]
  md$method <- ifelse(md$method == "CAMS", "Marginal-IPCW", "Marginal-AIPCW")
  md$R_label <- "all"
  diagnostics <- rbind(diagnostics, md)
  status <- data.frame(method = names(bounds), group = "all", status = "ok",
    message = support2_model_note(mdl, mdl0))
  if (include_subgroups) {
    bounds[["Separate-IPCW"]] <- NA_real_
    bounds[["Separate-AIPCW"]] <- NA_real_
    for (g in sort(unique(prepared$calib_groups))) {
      i <- prepared$fit_groups == g; j <- prepared$calib_groups == g
      k <- prepared$test_groups == g
      sub <- tryCatch({
        sx <- support2_subgroup_columns(fit[i, , drop = FALSE], xnames)
        if (!length(sx)) stop("No estimable subgroup predictors.")
        sm <- support2_checked_survreg(fit[i, , drop = FALSE], sx)
        sg <- support2_checked_survreg(fit[i, , drop = FALSE], sx, censoring = TRUE)
        if (any(!is.finite(coef(sm))) || any(!is.finite(coef(sg$fit)))) {
          stop("Nonfinite subgroup coefficients.")
        }
        calibrated <- est_alpha_ipcw(sm, test[k, sx, drop = FALSE], cal[j, , drop = FALSE], sx,
          alpha, sum(k), sg, probability_floor, bound_cap = cap,
          candidate_grid = candidate_grid, integration_method = "dense",
          integration_steps = integration_steps)
        calibrated$fit_note <- support2_model_note(sm, sg)
        calibrated
      }, error = function(e) e)
      if (inherits(sub, "error")) {
        status <- rbind(status, data.frame(method = c("Separate-IPCW", "Separate-AIPCW"),
          group = g, status = "failed", message = conditionMessage(sub)))
      } else {
        bounds[k, "Separate-IPCW"] <- sub$bounds$CAMS
        bounds[k, "Separate-AIPCW"] <- sub$bounds[["AIPCW-CAMS"]]
        sd <- sub$diagnostics[sub$diagnostics$method %in% c("CAMS", "AIPCW-CAMS"), ]
        sd$method <- ifelse(sd$method == "CAMS", "Separate-IPCW", "Separate-AIPCW")
        sd$R_label <- g
        diagnostics <- rbind(diagnostics, sd)
        status <- rbind(status, data.frame(method = c("Separate-IPCW", "Separate-AIPCW"),
          group = g, status = "ok", message = sub$fit_note))
      }
    }
  }
  if (include_cqr) {
    if (!requireNamespace("quantreg", quietly = TRUE)) {
      stop("Install quantreg, or run with --include-cqr=false.")
    }
    # Existing cqr() refers to rq() without a namespace.
    suppressPackageStartupMessages(library(quantreg))
  }
  baselines <- support2_baselines(prepared, mdl, mdl0, alpha, probability_floor,
    candidate_grid, include_subgroups, include_cqr, include_dft,
    dft_fixed_cutoff, imputation_seed)
  if (ncol(baselines$bounds)) bounds <- cbind(bounds, baselines$bounds)
  if (nrow(baselines$status)) status <- rbind(status, baselines$status)
  list(bounds = bounds, diagnostics = diagnostics, status = status,
    baselines = baselines, event_model = mdl, censoring_model = mdl0,
    runtime_seconds = unname(proc.time()[3] - start))
}

support2_evaluate <- function(prepared, result, probability_floor, integration_steps) {
  test <- prepared$data_test; ev <- prepared$evaluation
  cache <- make_aipcw_cache(result$event_model, data_calib = test,
    xnames = prepared$xnames, mdl0 = result$censoring_model,
    probability_floor = probability_floor, integration_method = "dense",
    bound_cap = prepared$horizon, integration_steps = integration_steps)
  groups <- c("all", sort(unique(ev$group)))
  rows <- list(); predictions <- list()
  for (method in names(result$bounds)) {
    L <- result$bounds[[method]]
    valid <- is.finite(L)
    phi <- rep(NA_real_, length(L)); loss <- phi; clipped <- phi
    if (any(valid)) {
      # Preserve scores for successful groups even if another subgroup fit failed.
      scored_cache <- if (all(valid)) cache else make_aipcw_cache(result$event_model,
        data_calib = test[valid, , drop = FALSE], xnames = prepared$xnames,
        mdl0 = result$censoring_model, probability_floor = probability_floor,
        integration_method = "dense", bound_cap = prepared$horizon,
        integration_steps = integration_steps)
      phi[valid] <- aipcw_risk_for_bound(L[valid], scored_cache, return_pseudo_loss = TRUE)
      loss[valid] <- test$event[valid] * as.numeric(test$censored_T[valid] < L[valid]) /
        scored_cache$G_y
      relevant_G <- as.numeric(sc_prob(result$censoring_model, test[valid, , drop = FALSE],
        prepared$xnames, pmin(test$censored_T[valid], L[valid])))
      clipped[valid] <- as.numeric(relevant_G < probability_floor)
    }
    for (g in groups) {
      idx <- if (g == "all") rep(TRUE, nrow(test)) else ev$group == g
      n <- sum(idx); complete <- all(valid[idx])
      exact_known <- ev$original_event[idx] == 1 | ev$original_Y[idx] >= L[idx]
      exact <- if (complete && all(exact_known)) mean(ev$original_Y[idx] >= L[idx]) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(method = method, audit_group = g,
        n_test = n, n_events = sum(test$event[idx]), complete = complete,
        coverage_exact = exact,
        coverage_ipcw = if (complete) 1 - mean(loss[idx]) else NA_real_,
        coverage_aipcw = if (complete) 1 - mean(phi[idx]) else NA_real_,
        coverage_observed_lower = if (complete) mean(test$censored_T[idx] >= L[idx]) else NA_real_,
        coverage_observed_upper = if (complete) 1 - mean(test$event[idx] == 1 &
          test$censored_T[idx] < L[idx]) else NA_real_,
        unknown_indicator_fraction = if (complete) mean(test$event[idx] == 0 &
          test$censored_T[idx] < L[idx]) else NA_real_,
        mean_bound = if (complete) mean(L[idx]) else NA_real_,
        median_bound = if (complete) median(L[idx]) else NA_real_,
        fraction_zero = if (complete) mean(L[idx] == 0) else NA_real_,
        fraction_at_cap = if (complete) mean(L[idx] == prepared$horizon) else NA_real_,
        relevant_G_clipped_fraction = if (complete) mean(clipped[idx]) else NA_real_)
    }
    predictions[[method]] <- data.frame(patient_id = ev$patient_id, group = ev$group,
      method = method, bound = L, original_Y = ev$original_Y,
      original_event = ev$original_event, observed_Y = test$censored_T,
      observed_event = test$event, pseudo_loss_ipcw = loss, pseudo_loss_aipcw = phi)
  }
  list(metrics = do.call(rbind, rows), predictions = do.call(rbind, predictions))
}

run_support2 <- function(repo_root = getwd(), csv_path = NULL, seed = 1L,
                         scenario = "native", horizon = 180, alpha = 0.1,
                         feature_set = "entry", landmark_days = 0,
                         probability_floor = 0.05, grid_step = 0.001,
                         integration_steps = 256L, include_subgroups = TRUE,
                         include_cqr = TRUE, output_dir = NULL,
                         include_dft = TRUE, dft_fixed_cutoff = "auto",
                         imputation_seed = NULL, aft_ridge_fallback = TRUE,
                         aft_ridge_theta = 0.1) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("Install the survival package.")
  suppressPackageStartupMessages(library(survival))
  source(file.path(repo_root, "utils/utils.R"), local = .GlobalEnv)
  source(file.path(repo_root, "utils/cams.R"), local = .GlobalEnv)
  source(file.path(repo_root, "utils/cqr.R"), local = .GlobalEnv)
  source(file.path(repo_root, "real_data/baselines.R"), local = .GlobalEnv)
  if (length(aft_ridge_fallback) != 1L || !is.logical(aft_ridge_fallback) ||
      is.na(aft_ridge_fallback) || !is.finite(aft_ridge_theta) || aft_ridge_theta <= 0)
    stop("Invalid AFT ridge policy.")
  old_options <- options(support2.aft_ridge_fallback = aft_ridge_fallback,
    support2.aft_ridge_theta = aft_ridge_theta)
  on.exit(options(old_options), add = TRUE)
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      !is.finite(probability_floor) || probability_floor <= 0 || probability_floor >= 1 ||
      !is.finite(grid_step) || grid_step <= 0 || grid_step > 0.1) stop("Invalid calibration setting.")
  if (is.null(imputation_seed)) imputation_seed <- support2_seed_offset(seed, 300001L)
  if (length(imputation_seed) != 1L || !is.finite(imputation_seed) ||
      imputation_seed < 0 || imputation_seed != as.integer(imputation_seed))
    stop("Invalid imputation seed.")
  if (length(integration_steps) != 1L || is.na(integration_steps) ||
      integration_steps < 16L) stop("integration_steps must be at least 16.")
  if (as.character(dft_fixed_cutoff) != "auto") {
    c0 <- suppressWarnings(as.numeric(dft_fixed_cutoff))
    if (!is.finite(c0) || c0 <= 0 || c0 > horizon) stop("Invalid DFT fixed cutoff.")
  }
  raw <- load_support2(csv_path, file.path(repo_root, "real_data/data"))
  prepared <- prepare_support2(raw, seed, scenario, horizon, feature_set, landmark_days)
  grid <- seq(grid_step, 1 - grid_step / 2, by = grid_step)
  grid <- grid[grid < 1]
  cat(sprintf("SUPPORT2: scenario=%s seed=%d, fit=%d cal=%d test=%d, predictors=%d\n",
    scenario, seed, nrow(prepared$data_fit), nrow(prepared$data_calib),
    nrow(prepared$data_test), length(prepared$xnames)))
  result <- support2_methods(prepared, alpha, probability_floor, grid,
    integration_steps, include_subgroups, include_cqr, include_dft,
    dft_fixed_cutoff, imputation_seed)
  evaluation <- support2_evaluate(prepared, result, probability_floor, integration_steps)
  evaluation$metrics$seed <- seed; evaluation$metrics$scenario <- scenario
  evaluation$metrics$target_coverage <- 1 - alpha
  evaluation$metrics$runtime_seconds <- result$runtime_seconds
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(evaluation$metrics, file.path(output_dir, "metrics.csv"), row.names = FALSE)
    write.csv(evaluation$predictions, file.path(output_dir, "predictions.csv"), row.names = FALSE)
    write.csv(result$diagnostics, file.path(output_dir, "calibration_diagnostics.csv"), row.names = FALSE)
    write.csv(result$status, file.path(output_dir, "method_status.csv"), row.names = FALSE)
    write.csv(result$baselines$diagnostics, file.path(output_dir, "baseline_diagnostics.csv"), row.names = FALSE)
    write.csv(result$baselines$imputation_manifest, file.path(output_dir, "dft_imputation_manifest.csv"), row.names = FALSE)
    method_names <- names(result$bounds)
    registry <- data.frame(method = method_names,
      scope = ifelse(grepl("subgroup|Separate", method_names), "subgroup fitting and calibration",
        ifelse(grepl("^CAMS", method_names), "pooled fitting; group calibration", "joint")),
      censoring_input = ifelse(grepl("DFT", method_names), "conditional-imputation adaptation",
        "original observed records"),
      candidate_learner = ifelse(grepl("CQR", method_names), "linear quantile regression of observed time",
        "Weibull AFT"),
      ct_learner = ifelse(grepl("DFT-CT", method_names), "lognormal AFT censoring quantile", ""))
    write.csv(registry, file.path(output_dir, "method_registry.csv"), row.names = FALSE)
    write.csv(prepared$split_manifest, file.path(output_dir, "split_manifest.csv"), row.names = FALSE)
    config <- list(seed = seed, scenario = scenario, horizon = horizon, alpha = alpha,
      feature_set = feature_set, landmark_days = landmark_days,
      probability_floor = probability_floor, grid_step = grid_step,
      integration_method = "dense", integration_steps = integration_steps,
      include_subgroups = include_subgroups, include_cqr = include_cqr,
      include_dft = include_dft, dft_fixed_cutoff = as.character(dft_fixed_cutoff),
      aft_ridge_fallback = aft_ridge_fallback, aft_ridge_theta = aft_ridge_theta,
      dft_censoring_input = "observed C for censoring; conditional lognormal imputation above Y for deaths",
      dft_ct_learner = "lognormal AFT conditional censoring quantile",
      imputation_seed = imputation_seed,
      code_md5 = as.list(tools::md5sum(file.path(repo_root, c("utils/cams.R", "utils/cqr.R",
        "utils/utils.R", "real_data/support2.R", "real_data/baselines.R")))),
      source_url = support2_source_url, source_md5 = attr(raw, "source_md5"),
      n_raw = nrow(raw), n_excluded_landmark = prepared$n_excluded_landmark,
      predictors = prepared$predictors, encoded_predictors = prepared$xnames,
      n_fit = nrow(prepared$data_fit), n_calib = nrow(prepared$data_calib), n_test = nrow(prepared$data_test),
      evaluation_note = "Exact coverage uses original observed outcomes; common fitted IPCW/AIPCW evaluators use fitting data only.")
    dput(config, file.path(output_dir, "config.R"))
    saveRDS(list(event_model = result$event_model, censoring_model = result$censoring_model,
      encoder = prepared$encoder, baseline_models = result$baselines$models,
      dft_cutoff_search = result$baselines$cutoff_search), file.path(output_dir, "models.rds"))
    capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  }
  list(metrics = evaluation$metrics, predictions = evaluation$predictions,
    diagnostics = result$diagnostics, status = result$status,
    prepared = prepared, result = result)
}

summarize_support2 <- function(metrics) {
  finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  combinations <- unique(metrics[, c("scenario", "method")])
  rows <- lapply(seq_len(nrow(combinations)), function(i) {
    scenario <- combinations$scenario[i]; method <- combinations$method[i]
    m <- metrics[metrics$scenario == scenario & metrics$method == method, , drop = FALSE]
    overall <- m[m$audit_group == "all", , drop = FALSE]
    successful <- overall$seed[overall$complete]
    overall_ok <- overall[overall$complete, , drop = FALSE]
    worst <- vapply(successful, function(seed) {
      g <- m[m$seed == seed & m$audit_group != "all", , drop = FALSE]
      if (nrow(g) && all(is.finite(g$coverage_exact))) min(g$coverage_exact) else NA_real_
    }, numeric(1))
    data.frame(scenario = scenario, method = method, n_runs = nrow(overall),
      n_complete_runs = nrow(overall_ok), n_exact_runs = sum(is.finite(overall_ok$coverage_exact)),
      mean_exact_coverage = finite_mean(overall_ok$coverage_exact),
      mean_worst_exact_group_coverage = finite_mean(worst),
      mean_ipcw_coverage = finite_mean(overall_ok$coverage_ipcw),
      mean_aipcw_coverage = finite_mean(overall_ok$coverage_aipcw),
      mean_bound = finite_mean(overall_ok$mean_bound))
  })
  do.call(rbind, rows)
}
