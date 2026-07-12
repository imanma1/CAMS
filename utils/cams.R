cams <- function(x, p, len_x, xnames,
                 data_fit, data_calib,
                 mdl0, alpha,
                 use_oracle_sc = FALSE) {

  newdata <- data.frame(x)
  colnames(newdata) <- xnames

  fmla <- as.formula(
    paste("Surv(censored_T, event) ~ ", paste(xnames, collapse = "+"))
  )

  mdl <- survreg(fmla, data = data_fit, dist = "weibull")

  eta <- 1 / log(nrow(data_calib))

  idx_test_0 <- newdata$X1 == 0
  idx_test_1 <- newdata$X1 == 1

  newdata0 <- newdata[idx_test_0, , drop = FALSE]
  newdata1 <- newdata[idx_test_1, , drop = FALSE]

  calib0 <- data_calib[data_calib$X1 == 0, , drop = FALSE]
  calib1 <- data_calib[data_calib$X1 == 1, , drop = FALSE]

  res0 <- est_alpha_ipcw(
    mdl, newdata0, calib0,
    xnames, alpha, nrow(newdata0), mdl0, eta, use_oracle_sc
  )

  res1 <- est_alpha_ipcw(
    mdl, newdata1, calib1,
    xnames, alpha, nrow(newdata1), mdl0, eta, use_oracle_sc
  )

  method_names <- names(res0)

  output <- as.data.frame(
    setNames(
      replicate(length(method_names), rep(NA_real_, nrow(newdata)), simplify = FALSE),
      method_names
    ),
    check.names = FALSE
  )

  for (method_name in method_names) {
    output[idx_test_0, method_name] <- res0[[method_name]]
    output[idx_test_1, method_name] <- res1[[method_name]]
  }

  output[] <- lapply(output, function(z) pmax(z, 0))

  return(list(
    output = output,
    times = rep(NA_real_, ncol(output))
  ))
}

est_alpha_ipcw <- function(mdl, newdata, data_calib,
                           xnames, alpha, len_x,
                           mdl0, eta,
                           use_oracle_sc = FALSE) {

  method_names <- c(
    "CAMS",
    "SN-CAMS",
    "CAMS-raw",
    "SN-CAMS-raw",
    "AIPCW-CAMS",
    "AIPCW-CAMS-raw"
  )

  if (nrow(newdata) == 0) {
    return(
      setNames(
        replicate(
          length(method_names),
          numeric(0),
          simplify = FALSE
        ),
        method_names
      )
    )
  }

  if (nrow(data_calib) == 0) {
    return(
      setNames(
        replicate(
          length(method_names),
          rep(0, nrow(newdata)),
          simplify = FALSE
        ),
        method_names
      )
    )
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
      return(
        rep(0, len_x)
      )
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

    lower_bnd
  }

  out <- lapply(
    method_names,
    get_bound_for_method
  )

  names(out) <- method_names

  out
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
      G_y = G_y,
      grid = numeric(0),
      cum_A = matrix(numeric(0), nrow = n, ncol = 0),
      cum_B = matrix(numeric(0), nrow = n, ncol = 0),
      mdl = mdl,
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

  S_T_grid <- survreg_weibull_survival(
    mdl = mdl,
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
    G_y = G_y,
    grid = censoring_grid,
    cum_A = row_cumsum(increment_A),
    cum_B = row_cumsum(increment_B),
    mdl = mdl,
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

    S_T_lower <- survreg_weibull_survival(
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