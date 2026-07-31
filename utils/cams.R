make_cams_R_label <- function(
    data,
    use_intersectional_R = FALSE
) {

  if (!"X1" %in% colnames(data)) {
    stop("X1 is required to define the CAMS groups.")
  }

  if (
    anyNA(data$X1) ||
      !all(data$X1 %in% c(0, 1))
  ) {
    stop("X1 must contain only 0 and 1.")
  }

  if (!use_intersectional_R) {
    return(
      ifelse(
        data$X1 == 0,
        "X1_0",
        "X1_1"
      )
    )
  }

  if (!"X2" %in% colnames(data)) {
    stop(
      "X2 is required when use_intersectional_R = TRUE."
    )
  }

  if (anyNA(data$X2)) {
    stop("X2 cannot contain missing values.")
  }

  paste0(
    "X1_", data$X1,
    "__X2_",
    ifelse(
      data$X2 > 0,
      "gt0",
      "le0"
    )
  )
}

make_oracle_event_model <- function(setting,
                                    homoscedastic_event = FALSE) {

  supported_settings <- c(
    "cov_cens_dr",
    "cov_cens_dr_tail",
    "cov_cens_dr_tail_hetero",
    "cams_vs_vanilla_lower_tail_hd_mild",
    "cams_vs_vanilla_lower_tail_hd_main",
    "cams_vs_vanilla_lower_tail_hd_strong"
  )

  if (!(setting %in% supported_settings)) {
    stop(
      sprintf(
        "Oracle event augmentation is not implemented for setting: %s",
        setting
      )
    )
  }

  structure(
    list(
      setting = setting,
      homoscedastic_event = homoscedastic_event
    ),
    class = "oracle_event_model"
  )
}


oracle_event_survival_prob <- function(mdl,
                                       newdata,
                                       t) {

  X <- as.data.frame(newdata)

  # ==================================================
  # Define the true event model for each setting
  # ==================================================

  if (
    mdl$setting %in% c(
      "cov_cens_dr",
      "cov_cens_dr_tail"
    )
  ) {

    # --------------------------------------------------
    # Simple constant-scale event models
    # --------------------------------------------------

    required_names <- c(
      "X1",
      "X2"
    )

    missing_names <- setdiff(
      required_names,
      colnames(X)
    )

    if (length(missing_names) > 0L) {
      stop(
        sprintf(
          "Oracle event model is missing columns: %s",
          paste(
            missing_names,
            collapse = ", "
          )
        )
      )
    }

    mu_t <- 2 +
      0.5 * X$X2 -
      0.5 * X$X1

    sigma_t <- rep(
      0.5,
      nrow(X)
    )


  } else if (
    mdl$setting == "cov_cens_dr_tail_hetero"
  ) {

    # --------------------------------------------------
    # Exact oracle event model for the heterogeneous
    # scale DR experiment.
    #
    # log(T) = mu_t(X) + sigma_t(X3) * epsilon_EV
    # --------------------------------------------------

    required_names <- c(
      "X1",
      "X2",
      "X3"
    )

    missing_names <- setdiff(
      required_names,
      colnames(X)
    )

    if (length(missing_names) > 0L) {
      stop(
        sprintf(
          "Oracle event model is missing columns: %s",
          paste(
            missing_names,
            collapse = ", "
          )
        )
      )
    }

    mu_t <- 2 +
      0.5 * X$X2 -
      0.5 * X$X1

    sigma_t <- ifelse(
      X$X3 > 0,
      0.60,
      0.35
    )


  } else if (
    mdl$setting %in% c(
      "cams_vs_vanilla_lower_tail_hd_mild",
      "cams_vs_vanilla_lower_tail_hd_main",
      "cams_vs_vanilla_lower_tail_hd_strong"
    )
  ) {

    # --------------------------------------------------
    # Existing lower-tail HD event model
    # --------------------------------------------------

    required_names <- paste0("X", 1:75)

    missing_names <- setdiff(
      required_names,
      colnames(X)
    )

    if (length(missing_names) > 0L) {
      stop(
        sprintf(
          "Oracle event model is missing columns: %s",
          paste(missing_names, collapse = ", ")
        )
      )
    }

    # Same dense signal used in the three lower-tail HD settings
    beta_dense <- 0.03 * rep(
      c(1, -1),
      length.out = 75 - 4
    )

    dense_score <- as.numeric(
      as.matrix(
        X[, paste0("X", 5:75), drop = FALSE]
      ) %*% beta_dense
    )

    mu_t <- 2.8 +
      0.40 * X$X1 +
      0.60 * X$X2 -
      0.50 * X$X3 +
      0.30 * X$X4 +
      dense_score

    if (isTRUE(mdl$homoscedastic_event)) {
      sigma_t <- rep(
        0.35,
        nrow(X)
      )
    } else {
      sigma_t <- 0.28 +
        0.17 * X$X1
    }

  } else {

    stop(
      sprintf(
        "Oracle event survival is not implemented for setting: %s",
        mdl$setting
      )
    )
  }

  # ==================================================
  # Compute true conditional event survival
  # ==================================================

  n <- nrow(X)

  # --------------------------------------------------
  # Matrix of times
  # --------------------------------------------------
  if (is.matrix(t)) {

    if (nrow(t) != n) {
      stop(
        sprintf(
          "For matrix t, nrow(t) must equal nrow(newdata): %d versus %d.",
          nrow(t),
          n
        )
      )
    }

    k <- ncol(t)

    mu_mat <- matrix(
      mu_t,
      nrow = n,
      ncol = k
    )

    sigma_mat <- matrix(
      sigma_t,
      nrow = n,
      ncol = k
    )

    t_safe <- pmax(
      t,
      .Machine$double.xmin
    )

    z <- (
      log(t_safe) -
        mu_mat
    ) / sigma_mat

    # Prevent exp(z) overflow
    z <- pmin(
      z,
      700
    )

    # For epsilon = log(-log(U)):
    # S_T(t | X) = exp(-exp(z))
    surv <- exp(
      -exp(z)
    )

    surv[t <= 0] <- 1

  } else {

    # --------------------------------------------------
    # Scalar or rowwise vector of times
    # --------------------------------------------------

    if (length(t) == 1L) {
      t <- rep(
        t,
        n
      )
    }

    if (length(t) != n) {
      stop(
        sprintf(
          "Length of t must equal nrow(newdata): %d versus %d.",
          length(t),
          n
        )
      )
    }

    t_safe <- pmax(
      t,
      .Machine$double.xmin
    )

    z <- (
      log(t_safe) -
        mu_t
    ) / sigma_t

    z <- pmin(
      z,
      700
    )

    surv <- exp(
      -exp(z)
    )

    surv[t <= 0] <- 1
  }

  surv[!is.finite(surv)] <- 0

  pmin(
    pmax(surv, 0),
    1
  )
}


augmentation_survival_prob <- function(mdl,
                                       newdata,
                                       t) {

  if (inherits(mdl, "oracle_event_model")) {
    return(
      oracle_event_survival_prob(
        mdl = mdl,
        newdata = newdata,
        t = t
      )
    )
  }

  # Existing fitted Weibull survreg augmentation model
  survreg_weibull_survival(
    mdl = mdl,
    newdata = newdata,
    t = t
  )
}

cams <- function(x, p, len_x, xnames,
                 data_fit, data_calib,
                 mdl0, alpha,
                 use_oracle_sc = FALSE,
                 augmentation_method = c(
                   "same",
                   "correct",
                   "wrong",
                   "oracle_event"
                 ),
                 augmentation_mdl = NULL,
                 setting = NULL,
                 homoscedastic_event = FALSE,
                 use_intersectional_R = FALSE) {

  newdata <- data.frame(x)
  colnames(newdata) <- xnames

  fmla <- as.formula(
    paste("Surv(censored_T, event) ~ ", paste(xnames, collapse = "+"))
  )

  mdl <- survreg(fmla, data = data_fit, dist = "weibull")

  augmentation_method <- match.arg(augmentation_method)

  if (is.null(augmentation_mdl)) {

    if (augmentation_method %in% c("same", "correct")) {

      # Practical full fitted Weibull model.
      #
      # "same" and the old "correct" retain their old behavior
      # so previous results remain interpretable.
      augmentation_mdl <- mdl

    } else if (augmentation_method == "oracle_event") {

      if (is.null(setting)) {
        stop(
          "setting must be supplied when augmentation_method = 'oracle_event'."
        )
      }

      augmentation_mdl <- make_oracle_event_model(
        setting = setting,
        homoscedastic_event = homoscedastic_event
      )

    } else if (augmentation_method == "wrong") {

      wrong_xnames <- intersect(
        c("X1", "X2"),
        xnames
      )

      if (length(wrong_xnames) == 0L) {
        stop(
          "The deliberately wrong augmentation model needs X1 or X2."
        )
      }

      wrong_fmla <- as.formula(
        paste(
          "Surv(censored_T, event) ~",
          paste(
            wrong_xnames,
            collapse = " + "
          )
        )
      )

      augmentation_mdl <- survreg(
        wrong_fmla,
        data = data_fit,
        dist = "weibull"
      )
    }
  }

  eta <- 1 / log(
    nrow(data_calib)
  )

  # ----------------------------------------------------------
  # Define the calibration groups R
  # ----------------------------------------------------------

  if (use_intersectional_R) {

    R_levels <- c(
      "X1_0__X2_le0",
      "X1_0__X2_gt0",
      "X1_1__X2_le0",
      "X1_1__X2_gt0"
    )

  } else {

    R_levels <- c(
      "X1_0",
      "X1_1"
    )
  }

  R_test <- make_cams_R_label(
    data = newdata,
    use_intersectional_R = use_intersectional_R
  )

  R_calib <- make_cams_R_label(
    data = data_calib,
    use_intersectional_R = use_intersectional_R
  )

  group_results <- setNames(
    vector(
      mode = "list",
      length = length(R_levels)
    ),
    R_levels
  )

  diagnostics_list <- setNames(
    vector(
      mode = "list",
      length = length(R_levels)
    ),
    R_levels
  )

  output <- NULL
  method_names <- NULL

  # ----------------------------------------------------------
  # Calibrate separately inside every R group
  # ----------------------------------------------------------

  for (r in R_levels) {

    idx_test_r <- R_test == r
    idx_calib_r <- R_calib == r

    n_test_r <- sum(idx_test_r)
    n_calib_r <- sum(idx_calib_r)

    cat(
      sprintf(
        paste0(
          "Calibrating CAMS for %s: ",
          "calibration=%d, test=%d\n"
        ),
        r,
        n_calib_r,
        n_test_r
      )
    )

    if (n_calib_r == 0L) {
      stop(
        sprintf(
          "No calibration observations are available for group %s.",
          r
        )
      )
    }

    if (n_test_r == 0L) {
      warning(
        sprintf(
          "No test observations are available for group %s.",
          r
        )
      )
    }

    res_r <- est_alpha_ipcw(
      mdl = mdl,
      newdata = newdata[
        idx_test_r,
        ,
        drop = FALSE
      ],
      data_calib = data_calib[
        idx_calib_r,
        ,
        drop = FALSE
      ],
      xnames = xnames,
      alpha = alpha,
      len_x = n_test_r,
      mdl0 = mdl0,
      eta = eta,
      use_oracle_sc = use_oracle_sc,
      augmentation_mdl = augmentation_mdl
    )

    group_results[[r]] <- res_r

    if (is.null(output)) {

      method_names <- names(
        res_r$bounds
      )

      output <- as.data.frame(
        setNames(
          replicate(
            length(method_names),
            rep(
              NA_real_,
              nrow(newdata)
            ),
            simplify = FALSE
          ),
          method_names
        ),
        check.names = FALSE
      )

    } else if (
      !identical(
        method_names,
        names(res_r$bounds)
      )
    ) {

      stop(
        sprintf(
          "CAMS method names differ for group %s.",
          r
        )
      )
    }

    for (method_name in method_names) {

      output[
        idx_test_r,
        method_name
      ] <- res_r$bounds[[
        method_name
      ]]
    }

    diagnostics_r <- res_r$diagnostics
    diagnostics_r$R_label <- r

    diagnostics_list[[r]] <- diagnostics_r
  }

  output[] <- lapply(
    output,
    function(z) {

      finite_index <- is.finite(z)

      z[finite_index] <- pmax(
        z[finite_index],
        0
      )

      z
    }
  )

  diagnostics <- do.call(
    rbind,
    diagnostics_list
  )

  rownames(diagnostics) <- NULL

  return(
    list(
      output = output,
      times = rep(
        NA_real_,
        ncol(output)
      ),
      diagnostics = diagnostics,
      augmentation_method = augmentation_method,
      R_levels = R_levels
    )
  )
}

est_alpha_ipcw <- function(mdl, newdata, data_calib,
                           xnames, alpha, len_x,
                           mdl0, eta,
                           use_oracle_sc = FALSE,
                           augmentation_mdl = mdl) {

  method_names <- c(
    "CAMS",
    "SN-CAMS",
    "CAMS-raw",
    "SN-CAMS-raw",
    "AIPCW-CAMS",
    "AIPCW-CAMS-raw"
  )

  if (nrow(newdata) == 0) {
    empty_bounds <- setNames(
        replicate(
          length(method_names),
          numeric(0),
          simplify = FALSE
        ),
        method_names
      )
    return(list(
      bounds = empty_bounds,
      diagnostics = data.frame(
        method = method_names,
        selected_calibration_level = NA_real_,
        estimated_calibration_risk = NA_real_,
        fraction_censoring_probabilities_truncated = NA_real_,
        n_calib = 0L,
        check.names = FALSE
      )
    ))
  }

  if (nrow(data_calib) == 0) {
    empty_bounds <- setNames(
        replicate(
          length(method_names),
          rep(0, nrow(newdata)),
          simplify = FALSE
        ),
        method_names
      )
    return(list(
      bounds = empty_bounds,
      diagnostics = data.frame(
        method = method_names,
        selected_calibration_level = 0,
        estimated_calibration_risk = 0,
        fraction_censoring_probabilities_truncated = NA_real_,
        n_calib = 0L,
        check.names = FALSE
      )
    ))
  }

  v_list <- seq(
    0.001,
    0.999,
    by = 0.001
  )

  calib_x <- data_calib[
    ,
    xnames,
    drop = FALSE
  ]

  n_calib_subgroup <- nrow(data_calib)

  # Stabilized AIPCW cache
  aipcw_cache <- make_aipcw_cache(
    mdl = mdl,
    augmentation_mdl = augmentation_mdl,
    data_calib = data_calib,
    xnames = xnames,
    mdl0 = mdl0,
    probability_floor = eta,
    use_oracle_sc = use_oracle_sc
  )

  # Approximately untruncated cache.
  # A tiny floor is retained only to prevent division by zero.
  aipcw_raw_cache <- make_aipcw_cache(
    mdl = mdl,
    augmentation_mdl = augmentation_mdl,
    data_calib = data_calib,
    xnames = xnames,
    mdl0 = mdl0,
    probability_floor = 1e-10,
    use_oracle_sc = use_oracle_sc
  )

  pr_calib <- aipcw_cache$G_y
  raw_pr_calib <- aipcw_raw_cache$G_y

  weight_calib <- 1 / pr_calib
  raw_weight_calib <- 1 / raw_pr_calib

  n_eff <- (
    n_calib_subgroup^2
  ) / sum(weight_calib^2)

  K_maps <- length(v_list)

  total_weight_event <- sum(
    weight_calib[
      data_calib$event == 1
    ]
  )

  total_raw_weight_event <- sum(
    raw_weight_calib[
      data_calib$event == 1
    ]
  )

  est_alpha_all <- function(v) {

    lv_calib <- as.numeric(
      lv_cams(
        mdl,
        calib_x,
        v
      )
    )

    ind <- (
      data_calib$censored_T < lv_calib
    ) & (
      data_calib$event == 1
    )

    sum_num <- sum(
      weight_calib[ind]
    )

    sum_raw_num <- sum(
      raw_weight_calib[ind]
    )

    # --------------------------------------------
    # IPCW estimators
    # --------------------------------------------

    risk_cams <- sum_num /
      n_calib_subgroup

    if (
      total_weight_event == 0 ||
      !is.finite(total_weight_event)
    ) {
      risk_sn_cams <- 1
    } else {
      risk_sn_cams <- sum_num /
        total_weight_event
    }

    risk_cams_raw <- sum_raw_num /
      n_calib_subgroup

    if (
      total_raw_weight_event == 0 ||
      !is.finite(total_raw_weight_event)
    ) {
      risk_sn_cams_raw <- 1
    } else {
      risk_sn_cams_raw <- sum_raw_num /
        total_raw_weight_event
    }

    # --------------------------------------------
    # AIPCW estimators
    # --------------------------------------------

    risk_aipcw <- aipcw_risk_for_bound(
      lower_bound = lv_calib,
      cache = aipcw_cache
    )

    risk_aipcw_raw <- aipcw_risk_for_bound(
      lower_bound = lv_calib,
      cache = aipcw_raw_cache
    )

    # --------------------------------------------
    # Optional finite-sample penalty
    # --------------------------------------------

    C_0 <- 0

    penalty <- C_0 * sqrt(
      (
        log(2) +
          log(K_maps)
      ) / n_eff
    )

    risk_vec <- c(
      "CAMS" = risk_cams + penalty,
      "SN-CAMS" = risk_sn_cams + penalty,
      "CAMS-raw" = risk_cams_raw + penalty,
      "SN-CAMS-raw" = risk_sn_cams_raw + penalty,
      "AIPCW-CAMS" = risk_aipcw + penalty,
      "AIPCW-CAMS-raw" = risk_aipcw_raw + penalty
    )

    risk_vec[
      !is.finite(risk_vec)
    ] <- 1

    risk_vec
  }

  risk_mat <- t(
    vapply(
      v_list,
      est_alpha_all,
      FUN.VALUE = setNames(
        numeric(length(method_names)),
        method_names
      )
    )
  )

  colnames(risk_mat) <- method_names

  # Enforce monotonic estimated risk as v becomes more aggressive
  for (method_name in method_names) {
    risk_mat[, method_name] <- cummax(
      risk_mat[, method_name]
    )
  }

  get_bound_for_method <- function(method_name) {

    feasible_idx <- which(
      risk_mat[, method_name] <= alpha
    )

    if (length(feasible_idx) == 0L) {
      return(list(
        bound = rep(0, len_x),
        selected_v = 0,
        selected_risk = 0
      ))
    }

    v_hat <- max(
      v_list[feasible_idx]
    )

    lower_bnd <- as.numeric(
      lv_cams(
        mdl,
        newdata,
        v_hat
      )
    )

    if (length(lower_bnd) == 1L) {
      lower_bnd <- rep(
        lower_bnd,
        len_x
      )
    }

    list(
      bound = lower_bnd,
      selected_v = v_hat,
      selected_risk = risk_mat[feasible_idx[length(feasible_idx)], method_name]
    )
  }

  selected <- lapply(
    method_names,
    get_bound_for_method
  )

  names(selected) <- method_names

  out <- lapply(selected, `[[`, "bound")

  fraction_truncated <- mean(
    !is.finite(aipcw_raw_cache$G_y_raw) |
      aipcw_raw_cache$G_y_raw < eta
  )

  diagnostics <- data.frame(
    method = method_names,
    selected_calibration_level = vapply(
      selected,
      `[[`,
      numeric(1),
      "selected_v"
    ),
    estimated_calibration_risk = vapply(
      selected,
      `[[`,
      numeric(1),
      "selected_risk"
    ),
    fraction_censoring_probabilities_truncated = rep(
      fraction_truncated,
      length(method_names)
    ),
    n_calib = rep(n_calib_subgroup, length(method_names)),
    check.names = FALSE
  )

  list(bounds = out, diagnostics = diagnostics)
}


lv_cams <- function(mdl, calib_x, v){
  if(length(v) == 0) {
    return(rep(0, times = nrow(calib_x)))
  }
  lv_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = v)
  return(lv_calib)
}

survreg_weibull_survival <- function(mdl, newdata, t) {

  if (!identical(mdl$dist, "weibull")) {
    stop("survreg_weibull_survival() currently supports dist = 'weibull' only.")
  }

  lp <- as.numeric(
    predict(mdl, newdata = newdata, type = "lp")
  )

  sigma <- as.numeric(mdl$scale)
  n <- length(lp)

  if (is.matrix(t)) {

    if (nrow(t) != n) {
      stop("For matrix t, nrow(t) must equal nrow(newdata).")
    }

    lp_mat <- matrix(
      lp,
      nrow = n,
      ncol = ncol(t)
    )

    t_safe <- pmax(t, .Machine$double.xmin)

    z <- (log(t_safe) - lp_mat) / sigma
    z <- pmin(z, 700)

    surv <- exp(-exp(z))
    surv[t <= 0] <- 1

  } else {

    if (length(t) == 1L) {
      t <- rep(t, n)
    }

    if (length(t) != n) {
      stop("Length of t must equal nrow(newdata).")
    }

    t_safe <- pmax(t, .Machine$double.xmin)

    z <- (log(t_safe) - lp) / sigma
    z <- pmin(z, 700)

    surv <- exp(-exp(z))
    surv[t <= 0] <- 1
  }

  surv[!is.finite(surv)] <- 0

  pmin(
    pmax(surv, 0),
    1
  )
}

row_cumsum <- function(mat) {

  mat <- as.matrix(mat)

  if (ncol(mat) == 0L) {
    return(mat)
  }

  if (ncol(mat) == 1L) {
    return(mat)
  }

  out <- t(
    apply(mat, 1, cumsum)
  )

  matrix(
    out,
    nrow = nrow(mat),
    ncol = ncol(mat)
  )
}

make_aipcw_cache <- function(mdl,
                             augmentation_mdl = mdl,
                             data_calib,
                             xnames,
                             mdl0,
                             probability_floor,
                             use_oracle_sc = FALSE) {

  Y <- as.numeric(data_calib$censored_T)
  event <- as.numeric(data_calib$event)
  calib_x <- data_calib[, xnames, drop = FALSE]

  n <- nrow(data_calib)

  G_y_raw <- sc_prob(
    mdl0 = mdl0,
    data = data_calib,
    xnames = xnames,
    t = Y
  )

  G_y <- pmax(
    as.numeric(G_y_raw),
    probability_floor
  )

  censoring_grid <- sort(
    unique(
      Y[
        event == 0 &
          is.finite(Y) &
          Y > 0
      ]
    )
  )

  Q <- length(censoring_grid)

  if (Q == 0L) {
    return(list(
      Y = Y,
      event = event,
      G_y_raw = as.numeric(G_y_raw),
      G_y = G_y,
      grid = numeric(0),
      cum_A = matrix(numeric(0), nrow = n, ncol = 0),
      cum_B = matrix(numeric(0), nrow = n, ncol = 0),
      mdl = augmentation_mdl,
      calib_x = calib_x
    ))
  }

  grid_mat <- matrix(
    censoring_grid,
    nrow = n,
    ncol = Q,
    byrow = TRUE
  )

  Y_mat <- matrix(
    Y,
    nrow = n,
    ncol = Q
  )

  min_time_mat <- matrix(
    pmin(Y_mat, grid_mat),
    nrow = n,
    ncol = Q
  )

  G_grid_raw <- sc_prob(
    mdl0 = mdl0,
    data = data_calib,
    xnames = xnames,
    t = grid_mat
  )

  G_min_raw <- sc_prob(
    mdl0 = mdl0,
    data = data_calib,
    xnames = xnames,
    t = min_time_mat
  )

  G_grid <- pmax(
    as.matrix(G_grid_raw),
    probability_floor
  )

  G_min <- pmax(
    as.matrix(G_min_raw),
    probability_floor
  )

  if (
    nrow(G_grid) != n ||
    ncol(G_grid) != Q ||
    nrow(G_min) != n ||
    ncol(G_min) != Q
  ) {
    stop("sc_prob() returned a matrix with unexpected dimensions.")
  }

  # Lambda_C(Y_i wedge u_k | X_i)
  Lambda_min <- -log(G_min)

  if (Q == 1L) {

    dLambda <- matrix(
      Lambda_min[, 1],
      nrow = n,
      ncol = 1
    )

  } else {

    dLambda <- cbind(
      Lambda_min[, 1],
      Lambda_min[, 2:Q, drop = FALSE] -
        Lambda_min[, 1:(Q - 1), drop = FALSE]
    )
  }

  # Numerical protection against tiny negative increments
  dLambda[dLambda < 0 & dLambda > -1e-10] <- 0

  dN <- matrix(
    0,
    nrow = n,
    ncol = Q
  )

  censored_idx <- which(event == 0)

  if (length(censored_idx) > 0L) {

    censoring_column <- match(
      Y[censored_idx],
      censoring_grid
    )

    valid <- !is.na(censoring_column)

    dN[
      cbind(
        censored_idx[valid],
        censoring_column[valid]
      )
    ] <- 1
  }

  # Estimated censoring martingale increments
  dM <- dN - dLambda

  S_T_grid <- augmentation_survival_prob(
    mdl = augmentation_mdl,
    newdata = calib_x,
    t = grid_mat
  )

  S_T_grid <- pmax(
    S_T_grid,
    1e-12
  )

  increment_A <- dM / G_grid

  increment_B <- dM / (
    G_grid * S_T_grid
  )

  list(
    Y = Y,
    event = event,
    G_y_raw = as.numeric(G_y_raw),
    G_y = G_y,
    grid = censoring_grid,
    cum_A = row_cumsum(increment_A),
    cum_B = row_cumsum(increment_B),
    mdl = augmentation_mdl,
    calib_x = calib_x
  )
}

aipcw_risk_for_bound <- function(lower_bound,
                                 cache) {

  lower_bound <- as.numeric(lower_bound)

  if (length(lower_bound) != length(cache$Y)) {
    stop("lower_bound and calibration data have different lengths.")
  }

  loss_observed <- as.numeric(
    cache$Y < lower_bound
  )

  ipcw_part <- (
    cache$event *
      loss_observed
  ) / cache$G_y

  if (length(cache$grid) == 0L) {

    pseudo_loss <- ipcw_part

  } else {

    # Number of censoring-grid points below each lower bound
    grid_index <- findInterval(
      lower_bound,
      cache$grid
    )

    A_value <- rep(
      0,
      length(lower_bound)
    )

    B_value <- rep(
      0,
      length(lower_bound)
    )

    has_grid_point <- grid_index > 0L

    if (any(has_grid_point)) {

      row_idx <- which(has_grid_point)

      A_value[has_grid_point] <- cache$cum_A[
        cbind(
          row_idx,
          grid_index[has_grid_point]
        )
      ]

      B_value[has_grid_point] <- cache$cum_B[
        cbind(
          row_idx,
          grid_index[has_grid_point]
        )
      ]
    }

    S_T_lower <- augmentation_survival_prob(
      mdl = cache$mdl,
      newdata = cache$calib_x,
      t = lower_bound
    )

    # Sum_k eta(u_k)/G(u_k) dM_C(u_k)
    augmentation <- A_value -
      S_T_lower * B_value

    pseudo_loss <- ipcw_part + augmentation
  }

  if (any(!is.finite(pseudo_loss))) {
    return(1)
  }

  mean(pseudo_loss)
}
