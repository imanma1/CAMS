new_cams <- function(x, p, len_x, xnames,
                     data_fit, data_calib,
                     mdl0, alpha,
                     use_oracle_sc = FALSE,
                     v_list = seq(0.001, 0.999, by = 0.001)) {

  newdata <- data.frame(x)
  colnames(newdata) <- xnames

  fmla <- as.formula(
    paste("Surv(censored_T, event) ~ ", paste(xnames, collapse = "+"))
  )

  mdl <- survreg(fmla, data = data_fit, dist = "weibull")

  eta <- 1 / log(nrow(data_calib))

  idx_test_0 <- newdata$X1 == 0
  idx_test_1 <- newdata$X1 == 1

  idx_calib_0 <- data_calib$X1 == 0
  idx_calib_1 <- data_calib$X1 == 1

  res0 <- new_est_alpha_ipcw_group(
    mdl = mdl,
    newdata = newdata[idx_test_0, , drop = FALSE],
    data_calib = data_calib[idx_calib_0, , drop = FALSE],
    xnames = xnames,
    alpha = alpha,
    mdl0 = mdl0,
    eta = eta,
    v_list = v_list,
    use_oracle_sc = use_oracle_sc
  )

  res1 <- new_est_alpha_ipcw_group(
    mdl = mdl,
    newdata = newdata[idx_test_1, , drop = FALSE],
    data_calib = data_calib[idx_calib_1, , drop = FALSE],
    xnames = xnames,
    alpha = alpha,
    mdl0 = mdl0,
    eta = eta,
    v_list = v_list,
    use_oracle_sc = use_oracle_sc
  )

  method_names <- colnames(res0$output)

  output <- as.data.frame(
    setNames(
      replicate(length(method_names), rep(NA_real_, nrow(newdata)), simplify = FALSE),
      method_names
    ),
    check.names = FALSE
  )

  for (method_name in method_names) {
    output[idx_test_0, method_name] <- res0$output[[method_name]]
    output[idx_test_1, method_name] <- res1$output[[method_name]]
  }

  output[] <- lapply(output, function(z) pmax(z, 0))

  return(list(
    output = output,
    times = rep(NA_real_, ncol(output)),
    selected_lambdas = list(
      "X1=0" = res0$selected_lambdas,
      "X1=1" = res1$selected_lambdas
    ),
    selected_risks = list(
      "X1=0" = res0$selected_risks,
      "X1=1" = res1$selected_risks
    ),
    diagnostics = list(
      "X1=0" = res0$cell_diagnostics,
      "X1=1" = res1$cell_diagnostics
    )
  ))
}

new_est_alpha_ipcw_group <- function(mdl, newdata, data_calib,
                                     xnames, alpha,
                                     mdl0, eta,
                                     v_list,
                                     use_oracle_sc = FALSE) {

  method_names <- c(
    "Local-CAMS",
    "Local-SN-CAMS",
    "Local-CAMS-raw",
    "Local-SN-CAMS-raw"
  )

  if (nrow(newdata) == 0) {
    output <- as.data.frame(
      setNames(
        replicate(length(method_names), numeric(0), simplify = FALSE),
        method_names
      ),
      check.names = FALSE
    )

    return(list(
      output = output,
      selected_lambdas = NULL,
      selected_risks = NULL
    ))
  }

  if (nrow(data_calib) == 0) {
    output <- as.data.frame(
      setNames(
        replicate(length(method_names), rep(0, nrow(newdata)), simplify = FALSE),
        method_names
      ),
      check.names = FALSE
    )

    return(list(
      output = output,
      selected_lambdas = NULL,
      selected_risks = NULL
    ))
  }

  calib_x <- data_calib[, xnames, drop = FALSE]

  if (use_oracle_sc || is_sc_model(mdl0)) {
    raw_pr_calib <- sc_prob(
      mdl0 = mdl0,
      data = data_calib,
      xnames = xnames,
      t = data_calib$censored_T
    )
  } else {
    calib_mat <- as.matrix(data_calib[, xnames, drop = FALSE])
    gpr_mean <- mdl0$predict(calib_mat)
    gpr_sd <- mdl0$predict(calib_mat, se.fit = TRUE)$se
    raw_pr_calib <- pnorm((-data_calib$censored_T - gpr_mean) / gpr_sd)
  }

  pr_calib <- pmax(raw_pr_calib, eta)
  weight_calib <- 1 / pr_calib

  raw_weight_calib <- 1 / raw_pr_calib

  lv_calib_mat <- predict(
    mdl,
    newdata = calib_x,
    type = "quantile",
    p = v_list
  )

  lv_calib_mat <- as.matrix(lv_calib_mat)

  if (nrow(lv_calib_mat) != nrow(data_calib) &&
      ncol(lv_calib_mat) == nrow(data_calib)) {
    lv_calib_mat <- t(lv_calib_mat)
  }

  lv_new_mat <- predict(
    mdl,
    newdata = newdata,
    type = "quantile",
    p = v_list
  )

  lv_new_mat <- as.matrix(lv_new_mat)

  if (nrow(lv_new_mat) != nrow(newdata) &&
      ncol(lv_new_mat) == nrow(newdata)) {
    lv_new_mat <- t(lv_new_mat)
  }

  idx_calib_le <- data_calib$X2 <= 0
  idx_calib_gt <- data_calib$X2 > 0

  idx_new_le <- newdata$X2 <= 0
  idx_new_gt <- newdata$X2 > 0

  get_cell_sums <- function(idx_cell, weights) {
    if (sum(idx_cell) == 0) {
      return(rep(0, length(v_list)))
    }

    lv_cell <- lv_calib_mat[idx_cell, , drop = FALSE]

    ind_cell <- (data_calib$censored_T[idx_cell] < lv_cell) &
      (data_calib$event[idx_cell] == 1)

    colSums(sweep(ind_cell, 1, weights[idx_cell], `*`))
  }

  get_cell_utility <- function(idx_cell) {
    if (sum(idx_cell) == 0) {
      return(rep(0, length(v_list)))
    }

    colSums(lv_calib_mat[idx_cell, , drop = FALSE])
  }

  num_le <- get_cell_sums(idx_calib_le, weight_calib)
  num_gt <- get_cell_sums(idx_calib_gt, weight_calib)

  raw_num_le <- get_cell_sums(idx_calib_le, raw_weight_calib)
  raw_num_gt <- get_cell_sums(idx_calib_gt, raw_weight_calib)

  util_le <- get_cell_utility(idx_calib_le)
  util_gt <- get_cell_utility(idx_calib_gt)

  n_R <- nrow(data_calib)

  total_weight_event <- sum(weight_calib[data_calib$event == 1])
  total_raw_weight_event <- sum(raw_weight_calib[data_calib$event == 1])

  num_mat <- outer(num_le, num_gt, "+")
  raw_num_mat <- outer(raw_num_le, raw_num_gt, "+")

  utility_mat <- outer(util_le, util_gt, "+") / n_R

  risk_cams <- num_mat / n_R
  risk_cams_raw <- raw_num_mat / n_R

  if (total_weight_event == 0 || is.na(total_weight_event)) {
    risk_sn_cams <- matrix(1, nrow = length(v_list), ncol = length(v_list))
  } else {
    risk_sn_cams <- num_mat / total_weight_event
  }

  if (total_raw_weight_event == 0 || is.na(total_raw_weight_event)) {
    risk_sn_cams_raw <- matrix(1, nrow = length(v_list), ncol = length(v_list))
  } else {
    risk_sn_cams_raw <- raw_num_mat / total_raw_weight_event
  }

  clean_risk <- function(risk_mat) {
    risk_mat[is.na(risk_mat) | is.nan(risk_mat) | !is.finite(risk_mat)] <- 1
    risk_mat
  }

  risk_list <- list(
    "Local-CAMS" = clean_risk(risk_cams),
    "Local-SN-CAMS" = clean_risk(risk_sn_cams),
    "Local-CAMS-raw" = clean_risk(risk_cams_raw),
    "Local-SN-CAMS-raw" = clean_risk(risk_sn_cams_raw)
  )

  safe_ratio <- function(num, den) {
    out <- num / den
    out[is.na(out) | is.nan(out) | !is.finite(out)] <- 1
    out
  }

  n_le <- sum(idx_calib_le)
  n_gt <- sum(idx_calib_gt)

  event_weight_le <- sum(weight_calib[idx_calib_le & data_calib$event == 1])
  event_weight_gt <- sum(weight_calib[idx_calib_gt & data_calib$event == 1])

  raw_event_weight_le <- sum(raw_weight_calib[idx_calib_le & data_calib$event == 1])
  raw_event_weight_gt <- sum(raw_weight_calib[idx_calib_gt & data_calib$event == 1])

  risk_le_cams <- safe_ratio(num_le, n_le)
  risk_gt_cams <- safe_ratio(num_gt, n_gt)

  risk_le_cams_raw <- safe_ratio(raw_num_le, n_le)
  risk_gt_cams_raw <- safe_ratio(raw_num_gt, n_gt)

  risk_le_sn_cams <- safe_ratio(num_le, event_weight_le)
  risk_gt_sn_cams <- safe_ratio(num_gt, event_weight_gt)

  risk_le_sn_cams_raw <- safe_ratio(raw_num_le, raw_event_weight_le)
  risk_gt_sn_cams_raw <- safe_ratio(raw_num_gt, raw_event_weight_gt)

  cell_feasible_list <- list(
    "Local-CAMS" = outer(
      risk_le_cams <= alpha,
      risk_gt_cams <= alpha,
      "&"
    ),

    "Local-SN-CAMS" = outer(
      risk_le_sn_cams <= alpha,
      risk_gt_sn_cams <= alpha,
      "&"
    ),

    "Local-CAMS-raw" = outer(
      risk_le_cams_raw <= alpha,
      risk_gt_cams_raw <= alpha,
      "&"
    ),

    "Local-SN-CAMS-raw" = outer(
      risk_le_sn_cams_raw <= alpha,
      risk_gt_sn_cams_raw <= alpha,
      "&"
    )
  )

  output <- as.data.frame(
    setNames(
      replicate(length(method_names), rep(0, nrow(newdata)), simplify = FALSE),
      method_names
    ),
    check.names = FALSE
  )

  selected_lambdas <- setNames(vector("list", length(method_names)), method_names)
  selected_risks <- setNames(vector("list", length(method_names)), method_names)

  for (method_name in method_names) {

    risk_mat <- risk_list[[method_name]]

    feasible_mat <- (risk_mat <= alpha) & cell_feasible_list[[method_name]]

    if (!any(feasible_mat, na.rm = TRUE)) {
      output[[method_name]] <- rep(0, nrow(newdata))
      selected_lambdas[[method_name]] <- NULL
      selected_risks[[method_name]] <- NULL
      next
    }

    utility_feasible <- utility_mat
    utility_feasible[!feasible_mat] <- -Inf

    best_idx <- which(utility_feasible == max(utility_feasible), arr.ind = TRUE)[1, ]

    k_le <- best_idx[1]
    k_gt <- best_idx[2]

    lower_new <- rep(NA_real_, nrow(newdata))

    lower_new[idx_new_le] <- lv_new_mat[idx_new_le, k_le]
    lower_new[idx_new_gt] <- lv_new_mat[idx_new_gt, k_gt]

    output[[method_name]] <- pmax(lower_new, 0)

    selected_lambdas[[method_name]] <- c(
      "X2<=0" = v_list[k_le],
      "X2>0" = v_list[k_gt]
    )

    selected_risks[[method_name]] <- list(
      "group_risk" = risk_mat[k_le, k_gt],
      "cell_risk_X2<=0" = c(
        "Local-CAMS" = risk_le_cams[k_le],
        "Local-SN-CAMS" = risk_le_sn_cams[k_le],
        "Local-CAMS-raw" = risk_le_cams_raw[k_le],
        "Local-SN-CAMS-raw" = risk_le_sn_cams_raw[k_le]
      )[method_name],
      "cell_risk_X2>0" = c(
        "Local-CAMS" = risk_gt_cams[k_gt],
        "Local-SN-CAMS" = risk_gt_sn_cams[k_gt],
        "Local-CAMS-raw" = risk_gt_cams_raw[k_gt],
        "Local-SN-CAMS-raw" = risk_gt_sn_cams_raw[k_gt]
      )[method_name]
    )
  }

  cell_diagnostics <- data.frame(
    cell = c("X2<=0", "X2>0"),
    n = c(sum(idx_calib_le), sum(idx_calib_gt)),
    events = c(
      sum(idx_calib_le & data_calib$event == 1),
      sum(idx_calib_gt & data_calib$event == 1)
    ),
    censoring_rate = c(
      mean(data_calib$event[idx_calib_le] == 0),
      mean(data_calib$event[idx_calib_gt] == 0)
    ),
    n_eff = c(
      sum(idx_calib_le)^2 / sum(weight_calib[idx_calib_le]^2),
      sum(idx_calib_gt)^2 / sum(weight_calib[idx_calib_gt]^2)
    ),
    max_weight = c(
      max(weight_calib[idx_calib_le]),
      max(weight_calib[idx_calib_gt])
    )
  )

  return(list(
    output = output,
    selected_lambdas = selected_lambdas,
    selected_risks = selected_risks,
    cell_diagnostics = cell_diagnostics
  ))
}
