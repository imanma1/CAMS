cams <- function(x, p, len_x, xnames,
                 data_fit, data_calib,
                 mdl0, alpha) {

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
    xnames, alpha, nrow(newdata0), mdl0, eta
  )

  res1 <- est_alpha_ipcw(
    mdl, newdata1, calib1,
    xnames, alpha, nrow(newdata1), mdl0, eta
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
                           mdl0, eta) {

  method_names <- c(
    "CAMS",
    "SN-CAMS",
    "CAMS-raw",
    "SN-CAMS-raw"
  )

  if (nrow(newdata) == 0) {
    return(setNames(
      replicate(length(method_names), numeric(0), simplify = FALSE),
      method_names
    ))
  }

  if (nrow(data_calib) == 0) {
    return(setNames(
      replicate(length(method_names), rep(0, nrow(newdata)), simplify = FALSE),
      method_names
    ))
  }

  v_list <- seq(0.001, 0.999, by = 0.001)

  calib_mat <- as.matrix(data_calib[, xnames, drop = FALSE])

  gpr_mean <- mdl0$predict(calib_mat)
  gpr_sd <- mdl0$predict(calib_mat, se.fit = TRUE)$se

  calib_x <- data_calib[, xnames, drop = FALSE]
  n_calib_subgroup <- nrow(data_calib)

  raw_pr_calib <- pnorm((-data_calib$censored_T - gpr_mean) / gpr_sd)

  # Eta-truncated probabilities and weights
  pr_calib <- pmax(raw_pr_calib, eta)
  weight_calib <- 1 / pr_calib

  # Raw probabilities and raw weights: no eta truncation
  raw_weight_calib <- 1 / raw_pr_calib

  n_eff <- (n_calib_subgroup^2) / sum(weight_calib^2)
  K_maps <- length(v_list)

  total_weight_event <- sum(weight_calib[data_calib$event == 1])
  total_raw_weight_event <- sum(raw_weight_calib[data_calib$event == 1])

  est_alpha_all <- function(v) {
    lv_calib <- lv_cams(mdl, calib_x, v)

    ind <- (data_calib$censored_T < lv_calib) &
      (data_calib$event == 1)

    # Numerators
    sum_num <- sum(weight_calib[ind])
    sum_raw_num <- sum(raw_weight_calib[ind])

    # 1. Regular CAMS: eta-truncated weights, n_R denominator
    risk_cams <- sum_num / n_calib_subgroup

    # 2. SN-CAMS: eta-truncated weights, self-normalized denominator
    if (total_weight_event == 0 || is.na(total_weight_event)) {
      risk_sn_cams <- 1
    } else {
      risk_sn_cams <- sum_num / total_weight_event
    }

    # 3. CAMS-raw: raw weights, n_R denominator
    risk_cams_raw <- sum_raw_num / n_calib_subgroup

    # 4. SN-CAMS-raw: raw weights, self-normalized denominator
    if (total_raw_weight_event == 0 || is.na(total_raw_weight_event)) {
      risk_sn_cams_raw <- 1
    } else {
      risk_sn_cams_raw <- sum_raw_num / total_raw_weight_event
    }

    C_0 <- 0
    penalty <- C_0 * sqrt((log(2) + log(K_maps)) / n_eff)

    risk_cams <- risk_cams + penalty
    risk_sn_cams <- risk_sn_cams + penalty
    risk_cams_raw <- risk_cams_raw + penalty
    risk_sn_cams_raw <- risk_sn_cams_raw + penalty

    if (is.na(risk_cams) || is.nan(risk_cams)) risk_cams <- 1
    if (is.na(risk_sn_cams) || is.nan(risk_sn_cams)) risk_sn_cams <- 1
    if (is.na(risk_cams_raw) || is.nan(risk_cams_raw)) risk_cams_raw <- 1
    if (is.na(risk_sn_cams_raw) || is.nan(risk_sn_cams_raw)) risk_sn_cams_raw <- 1

    c(
      "CAMS" = risk_cams,
      "SN-CAMS" = risk_sn_cams,
      "CAMS-raw" = risk_cams_raw,
      "SN-CAMS-raw" = risk_sn_cams_raw
    )
  }

  risk_mat <- t(sapply(v_list, est_alpha_all))

  for (method_name in method_names) {
    risk_mat[, method_name] <- cummax(risk_mat[, method_name])
  }

  get_bound_for_method <- function(method_name) {
    feasible <- v_list[risk_mat[, method_name] <= alpha]

    if (length(feasible) == 0) {
      return(rep(0, len_x))
    }

    v_hat <- max(feasible)

    lower_bnd <- as.numeric(lv_cams(mdl, newdata, v_hat))

    if (length(lower_bnd) == 1) {
      lower_bnd <- rep(lower_bnd, len_x)
    }

    lower_bnd
  }

  out <- lapply(method_names, get_bound_for_method)
  names(out) <- method_names

  return(out)
}


lv_cams <- function(mdl, calib_x, v){
  if(length(v) == 0) {
    return(rep(0, times = nrow(calib_x)))
  }
  lv_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = v)
  return(lv_calib)
}