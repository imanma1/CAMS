# DFT formulas from the simulation code, adapted to general right censoring.
# Native outcomes are never replaced; only missing C is imputed for DFT.
# No package installation, forking, or inspection of test outcomes occurs here.

support2_capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = value, warnings = unique(warnings))
}

support2_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  force(expr)
}

support2_seed_offset <- function(seed, offset) {
  as.integer((as.double(seed) + offset) %% (.Machine$integer.max - 1))
}

support2_checked_survreg <- function(data, xnames, censoring = FALSE,
    ridge_fallback = getOption("support2.aft_ridge_fallback", TRUE),
    ridge_theta = getOption("support2.aft_ridge_theta", 0.1)) {
  dat <- data
  if (censoring) dat$event <- 1L - dat$event
  if (!any(dat$event == 1)) stop("No observed target events for this working model.")
  rhs <- if (length(xnames)) paste(xnames, collapse = " + ") else "1"
  form <- as.formula(paste("survival::Surv(censored_T, event) ~", rhs), env = .GlobalEnv)
  attempt <- function(f) support2_capture(survival::survreg(f, dat,
    dist = if (censoring) "lognormal" else "weibull",
    control = survival::survreg.control(maxiter = 1000)))
  usable <- function(z) !inherits(z$value, "error") &&
    !any(grepl("converg|iterations", z$warnings, ignore.case = TRUE)) &&
    all(is.finite(coef(z$value))) && length(z$value$scale) == 1L &&
    is.finite(z$value$scale) && z$value$scale > 0
  z <- attempt(form); note <- ""
  if (!usable(z)) {
    cause <- if (inherits(z$value, "error")) conditionMessage(z$value)
      else if (any(!is.finite(coef(z$value)))) paste("nonfinite coefficients:",
        paste(names(coef(z$value))[!is.finite(coef(z$value))], collapse = ","))
      else paste(z$warnings, collapse = " | ")
    if (!ridge_fallback || !length(xnames)) stop("Unusable unpenalized AFT fit: ", cause)
    ridge_form <- as.formula(paste("survival::Surv(censored_T, event) ~ survival::ridge(",
      paste(xnames, collapse = ","), ", theta =", ridge_theta, ")"), env = .GlobalEnv)
    z <- attempt(ridge_form)
    if (!usable(z)) stop("Unusable AFT fit even with ridge fallback: ", cause)
    note <- paste0("training-only ridge fallback theta=", ridge_theta, "; ", cause)
  }
  mdl <- z$value
  attr(mdl, "support2_fit_note") <- note
  out <- if (censoring) structure(list(fit = mdl, xnames = xnames),
    class = "sc_aft_lognormal") else mdl
  attr(out, "support2_fit_note") <- note
  out
}

support2_model_note <- function(...) {
  notes <- vapply(list(...), function(m) {
    note <- attr(m, "support2_fit_note"); if (is.null(note)) "" else note
  }, character(1))
  paste(unique(notes[nzchar(notes)]), collapse = " | ")
}

support2_impute_c <- function(data, censoring_model, seed) {
  if (!inherits(censoring_model, "sc_aft_lognormal"))
    stop("Conditional imputation currently requires lognormal AFT censoring.")
  Y <- data$censored_T; D <- data$event
  if (any(!is.finite(Y) | Y <= 0) || any(!D %in% c(0, 1)))
    stop("Invalid observed survival records for imputation.")
  C <- Y
  idx <- which(D == 1L)
  if (length(idx)) {
    mu <- as.numeric(predict(censoring_model$fit, data[idx, , drop = FALSE], type = "lp"))
    sigma <- censoring_model$fit$scale
    log_G_y <- pnorm((log(Y[idx]) - mu) / sigma, lower.tail = FALSE, log.p = TRUE)
    u <- support2_with_seed(seed, runif(length(idx)))
    # Log survival inversion remains stable when F_C(Y) rounds to one.
    z <- qnorm(log_G_y + log(pmax(u, .Machine$double.xmin)),
      lower.tail = FALSE, log.p = TRUE)
    C[idx] <- exp(mu + sigma * z)
  }
  if (any(!is.finite(C)) || any(C < Y) || any(C[D == 1L] <= Y[D == 1L]))
    stop("Conditional censoring imputation produced invalid times.")
  C
}

support2_dft_clip <- function(L, cap, n) {
  if (length(L) != n || anyNA(L) || any(L == Inf)) stop("Invalid baseline predictions.")
  # A conformal correction of Inf gives -Inf, a legitimate zero lower bound.
  pmin(cap, pmax(0, as.numeric(L)))
}

support2_dft_adaptive <- function(event_model, censoring_model, cal, test,
                                 C_cal, xnames, alpha, probability_floor,
                                 candidate_grid, cap, ct = FALSE) {
  if (length(C_cal) != nrow(cal) || any(!is.finite(C_cal)) ||
      any(C_cal < cal$censored_T)) stop("DFT requires coherent censoring times.")
  quantiles <- function(data, levels) {
    q <- as.matrix(predict(event_model, data, type = "quantile", p = levels))
    if (nrow(q) != nrow(data)) q <- t(q)
    if (nrow(q) != nrow(data) || ncol(q) != length(levels) || any(!is.finite(q)))
      stop("Invalid event quantile matrix.")
    q
  }
  censor_cap <- function(data) {
    if (!ct) return(rep(cap, nrow(data)))
    mu <- as.numeric(predict(censoring_model$fit, data, type = "lp"))
    pmin(cap, exp(mu + censoring_model$fit$scale * qnorm(1 - probability_floor)))
  }
  L <- pmin(quantiles(cal, candidate_grid), censor_cap(cal))
  G <- sc_prob(censoring_model, cal, xnames, L)
  W <- 1 / pmax(G, probability_floor)
  eligible <- C_cal >= L
  denominator <- colSums(W * eligible)
  risk <- colSums(W * eligible * (cal$censored_T < L)) / denominator
  risk[!is.finite(risk) | denominator <= 0] <- Inf
  # Increasing lambda corresponds to decreasing v in the historical code.
  # This is its monot() envelope with the candidate order reversed.
  envelope <- cummax(risk)
  feasible <- which(envelope <= alpha)
  if (!length(feasible)) return(list(bound = rep(0, nrow(test)),
    level = 0, risk = 0, n_selected = nrow(cal), effective_n = nrow(cal),
    clipped_fraction = 0, fallback = "no feasible candidate"))
  j <- max(feasible); level <- candidate_grid[j]
  w <- W[eligible[, j], j]
  list(bound = support2_dft_clip(pmin(as.numeric(quantiles(test, level)),
      censor_cap(test)), cap, nrow(test)), level = level, risk = risk[j],
    n_selected = length(w), effective_n = sum(w)^2 / sum(w^2),
    clipped_fraction = mean(G[, j] < probability_floor), fallback = "")
}

support2_dft_fixed <- function(event_model, censoring_model, cal, test,
                              C_cal, xnames, alpha, probability_floor,
                              cutoff, cap) {
  if (length(C_cal) != nrow(cal) || any(!is.finite(C_cal)) ||
      any(C_cal < cal$censored_T)) stop("DFT requires coherent censoring times.")
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 || cutoff > cap)
    stop("DFT fixed cutoff must be positive and no greater than the horizon.")
  keep <- C_cal >= cutoff
  if (!any(keep)) return(list(bound = rep(0, nrow(test)), cutoff = cutoff,
    correction_mean = Inf, n_selected = 0, effective_n = 0,
    clipped_fraction = NA_real_, fallback = "empty cutoff-selected calibration set"))
  selected <- cal[keep, , drop = FALSE]
  G_cal <- as.numeric(sc_prob(censoring_model, selected, xnames, cutoff))
  G_test <- as.numeric(sc_prob(censoring_model, test, xnames, cutoff))
  weights <- 1 / pmax(G_cal, probability_floor)
  weights_test <- 1 / pmax(G_test, probability_floor)
  q_cal <- as.numeric(predict(event_model, selected, type = "quantile", p = alpha))
  q_test <- as.numeric(predict(event_model, test, type = "quantile", p = alpha))
  if (any(!is.finite(c(q_cal, q_test)))) stop("Invalid fixed-DFT quantiles.")
  score <- pmin(cutoff, q_cal) - pmin(cutoff, selected$censored_T)
  ordering <- order(score)
  cumulative <- cumsum(weights[ordering])
  target <- (1 - alpha) * (sum(weights) + weights_test)
  # Test-point weight is placed at +Inf, as in get_calibration().
  index <- findInterval(target, cumulative, left.open = TRUE) + 1L
  correction <- c(score[ordering], Inf)[index]
  L <- support2_dft_clip(pmin(cutoff, q_test) - correction, cap, nrow(test))
  list(bound = pmin(L, cutoff), cutoff = cutoff, correction_mean = mean(correction),
    n_selected = sum(keep), effective_n = sum(weights)^2 / sum(weights^2),
    clipped_fraction = mean(c(G_cal, G_test) < probability_floor), fallback = "")
}

support2_select_dft_cutoff <- function(fit, C_fit, xnames, censoring_model,
                                      alpha, probability_floor, cap, seed,
                                      requested = "auto") {
  if (!identical(as.character(requested), "auto")) {
    c0 <- suppressWarnings(as.numeric(requested))
    if (length(c0) != 1L || !is.finite(c0) || c0 <= 0 || c0 > cap)
      stop("--dft-fixed-cutoff must be auto or a positive number <= horizon.")
    return(list(cutoff = c0, candidates = data.frame(cutoff = c0, mean_bound = NA_real_)))
  }
  # Tuning is confined to the outer fitting split; no outer calibration/test labels.
  n <- nrow(fit)
  indices <- support2_with_seed(seed, sample.int(n))
  i <- indices[seq_len(floor(n / 2))]
  j <- indices[seq.int(length(i) + 1L, length(i) + floor(n / 4))]
  k <- setdiff(indices, c(i, j))
  sx <- support2_subgroup_columns(fit[i, , drop = FALSE], xnames)
  inner_model <- support2_checked_survreg(fit[i, , drop = FALSE], sx)
  candidates <- sort(unique(pmin(cap, c(cap, as.numeric(quantile(C_fit,
    probs = seq(0.1, 0.9, by = 0.05), names = FALSE, type = 1))))))
  means <- vapply(candidates, function(c0) {
    z <- support2_capture(support2_dft_fixed(inner_model, censoring_model,
      fit[j, , drop = FALSE], fit[k, , drop = FALSE], C_fit[j], sx,
      alpha, probability_floor, c0, cap))$value
    if (inherits(z, "error")) -Inf else mean(z$bound)
  }, numeric(1))
  valid <- which(is.finite(means))
  if (!length(valid)) stop("All training-only DFT cutoff candidates failed.")
  best <- valid[which.max(means[valid])]
  list(cutoff = candidates[best], candidates = data.frame(cutoff = candidates, mean_bound = means))
}

support2_baselines <- function(prepared, event_model, censoring_model, alpha,
                               probability_floor, candidate_grid,
                               include_subgroups = TRUE, include_cqr = TRUE,
                               include_dft = TRUE, fixed_cutoff = "auto",
                               imputation_seed = 300001L) {
  n <- nrow(prepared$data_test)
  bounds <- data.frame(row.names = seq_len(n)); statuses <- list(); diagnostics <- list()
  manifests <- list(); models <- list(); cutoffs <- list()
  scopes <- c("joint", if (include_subgroups) "subgroup")
  for (scope in scopes) {
    groups <- if (scope == "joint") "all" else sort(unique(prepared$fit_groups))
    base_names <- c(if (include_dft) c("Adaptive DFT-T", "Adaptive DFT-CT", "DFT-fixed"),
      if (include_cqr) "CQR")
    for (name in base_names) bounds[[paste0(name, " (", scope, ")")]] <- NA_real_
    for (gg in seq_along(groups)) {
      g <- groups[gg]
      idx <- function(labels) if (scope == "joint") rep(TRUE, length(labels)) else labels == g
      i <- idx(prepared$fit_groups); j <- idx(prepared$calib_groups); k <- idx(prepared$test_groups)
      fit <- prepared$data_fit[i, , drop = FALSE]; cal <- prepared$data_calib[j, , drop = FALSE]
      test <- prepared$data_test[k, , drop = FALSE]
      sx <- if (scope == "joint") prepared$xnames else support2_subgroup_columns(fit, prepared$xnames)
      model_note <- ""
      run_method <- function(name, expr) {
        method <- paste0(name, " (", scope, ")"); started <- proc.time()[3]
        captured <- support2_capture(expr); z <- captured$value
        if (!inherits(z, "error")) {
          checked <- support2_capture(support2_dft_clip(z$bound, prepared$horizon, sum(k)))
          if (inherits(checked$value, "error")) z <- checked$value else z$bound <- checked$value
        }
        bad <- inherits(z, "error")
        if (!bad) bounds[k, method] <<- z$bound
        notes <- c(if (bad) conditionMessage(z), captured$warnings,
          if (!bad) z$fallback, model_note)
        note <- paste(notes[nzchar(notes)], collapse = " | ")
        statuses[[length(statuses) + 1L]] <<- data.frame(method = method, group = g,
          status = if (bad) "failed" else "ok", message = note)
        value <- function(field) if (bad || is.null(z[[field]])) NA_real_ else z[[field]]
        diagnostics[[length(diagnostics) + 1L]] <<- data.frame(method = method, group = g,
          scope = scope, status = if (bad) "failed" else "ok", selected_level = value("level"),
          selected_cutoff = value("cutoff"), estimated_risk = value("risk"),
          n_selected = value("n_selected"), effective_n = value("effective_n"),
          clipped_fraction = value("clipped_fraction"), correction_mean = value("correction_mean"),
          elapsed_seconds = unname(proc.time()[3] - started), message = note)
      }
      if (include_cqr) run_method("CQR", {
        combined <- rbind(fit, cal)
        qr <- cqr(test[, sx, drop = FALSE], combined[, sx, drop = FALSE],
          combined$censored_T, alpha, I_fit = seq_len(nrow(fit)), exact_quantile = TRUE)
        list(bound = qr, n_selected = nrow(cal), effective_n = nrow(cal), fallback = "")
      })
      if (!include_dft) next
      imputed <- support2_capture({
        sm <- if (scope == "joint") event_model else support2_checked_survreg(fit, sx)
        sg <- if (scope == "joint") censoring_model else support2_checked_survreg(fit, sx, TRUE)
        rng_seed <- support2_seed_offset(imputation_seed, if (scope == "joint") 0 else 10000L + gg * 100L)
        Cf <- support2_impute_c(fit, sg, rng_seed)
        Cc <- support2_impute_c(cal, sg, support2_seed_offset(rng_seed, 1L))
        list(event_model = sm, censoring_model = sg, C_fit = Cf, C_cal = Cc)
      })
      z <- imputed$value
      if (inherits(z, "error")) {
        for (name in c("Adaptive DFT-T", "Adaptive DFT-CT", "DFT-fixed"))
          run_method(name, stop(conditionMessage(z)))
        next
      }
      key <- paste(scope, g, sep = "__"); models[[key]] <- z[c("event_model", "censoring_model")]
      model_note <- support2_model_note(z$event_model, z$censoring_model)
      m <- prepared$split_manifest
      for (part in c("fit", "calibration")) {
        ii <- m$split == part & (scope == "joint" | m$group == g)
        manifests[[length(manifests) + 1L]] <- data.frame(patient_id = m$patient_id[ii],
          split = part, group = m$group[ii], scope = scope, original_Y = m$observed_Y[ii],
          original_event = m$observed_event[ii], C_used = if (part == "fit") z$C_fit else z$C_cal,
          is_imputed = m$observed_event[ii] == 1L)
      }
      for (ct in c(FALSE, TRUE)) run_method(if (ct) "Adaptive DFT-CT" else "Adaptive DFT-T",
        support2_dft_adaptive(z$event_model, z$censoring_model, cal, test, z$C_cal,
          sx, alpha, probability_floor, candidate_grid, prepared$horizon, ct))
      run_method("DFT-fixed", {
        selected <- support2_select_dft_cutoff(fit, z$C_fit, sx, z$censoring_model,
          alpha, probability_floor, prepared$horizon,
          support2_seed_offset(imputation_seed, 20000L + gg), fixed_cutoff)
        cutoffs[[key]] <- selected$candidates
        support2_dft_fixed(z$event_model, z$censoring_model, cal, test, z$C_cal, sx,
          alpha, probability_floor, selected$cutoff, prepared$horizon)
      })
    }
  }
  bind <- function(rows) if (length(rows)) do.call(rbind, rows) else data.frame()
  list(bounds = bounds, status = bind(statuses), diagnostics = bind(diagnostics),
    imputation_manifest = bind(manifests), models = models, cutoff_search = cutoffs)
}
