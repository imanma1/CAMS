############################################
## estimate miscoverage rate
## using estimated quantile of T
############################################
alpha_qt <- function(mdl, newdata, data_fit, data_calib, xnames, alpha, len_x, mdl0, cens_rt,
                     use_oracle_sc = FALSE) {
  v_list <- v_pts_qt(mdl, data_fit, data_calib, xnames, alpha, cens_rt)
  v_list <- sort(unique(as.numeric(v_list)))

  calib_x <- data_calib[, xnames, drop = FALSE]

  # === VECTORIZED THRESHOLD EVALUATION ===
  lv_mat <- predict(mdl,
                    newdata = calib_x,
                    type = "quantile",
                    p = 1 - v_list)

  lv_mat <- as.matrix(lv_mat)

  if (nrow(lv_mat) != nrow(data_calib) && ncol(lv_mat) == nrow(data_calib)) {
    lv_mat <- t(lv_mat)
  }

  if (use_oracle_sc || is_sc_model(mdl0)) {
    pr_mat <- sc_prob(
      mdl0 = mdl0,
      data = data_calib,
      xnames = xnames,
      t = lv_mat
    )
  } else {
    calib_mat <- as.matrix(data_calib[, xnames, drop = FALSE])
    gpr_mean <- mdl0$predict(calib_mat)
    gpr_sd <- mdl0$predict(calib_mat, se.fit = TRUE)$se
    pr_mat <- pnorm((-lv_mat - gpr_mean) / gpr_sd)
  }
  weight_mat <- 1 / pr_mat

  ind1_mat <- (data_calib$censored_T < lv_mat) & (data_calib$C >= lv_mat)
  ind2_mat <- data_calib$C >= lv_mat

  sum_num <- colSums(weight_mat * ind1_mat)
  sum_den <- colSums(weight_mat * ind2_mat)

  alpha_v_list <- sum_num / sum_den

  bad_alpha <- is.na(alpha_v_list) |
    is.nan(alpha_v_list) |
    !is.finite(alpha_v_list) |
    is.na(sum_den) |
    is.nan(sum_den) |
    !is.finite(sum_den) |
    sum_den == 0

  alpha_v_list[bad_alpha] <- 1

  # monotonize alpha
  alpha_v <- monot(alpha_v_list)

  bad_alpha_v <- is.na(alpha_v) |
    is.nan(alpha_v) |
    !is.finite(alpha_v)

  alpha_v[bad_alpha_v] <- 1

  feasible_idx <- which(alpha_v <= alpha)

  if (length(feasible_idx) == 0) {
    v_hat_l <- NULL
  } else {
    v_hat_l <- min(v_list[feasible_idx])
  }

  # === OPTIMIZATION 2: VECTORIZE PREDICTION ===
  if (is.null(v_hat_l)) {
    lower_bnd_l <- rep(0, len_x)
  } else {
    lower_bnd_l <- as.numeric(lv_qt(mdl, newdata, v_hat_l, alpha, cens_rt))

    # Safety net: If the function somehow still returns a scalar, expand it
    if (length(lower_bnd_l) == 1) {
      lower_bnd_l <- rep(lower_bnd_l, len_x)
    }
  }
  lower_bnd_g <- rep(0, len_x)

  return(list(lower_bnd_l = lower_bnd_l, lower_bnd_g = lower_bnd_g))
}


################################################################
## The function returns predictive intervals resulting
## from the L_v defined based on integrated quantiles
################################################################
alpha_qct <- function(mdl, qc_mdl, newdata, data_fit, data_calib, xnames, alpha, len_x, mdl0, cens_rt,
                      use_oracle_sc = FALSE){
  v_list <- v_pts_qct(mdl, qc_mdl, data_fit, data_calib, xnames, alpha, cens_rt)
  v_list <- sort(unique(as.numeric(v_list)))

  calib_x <- data_calib[, xnames, drop = FALSE]

  # === OPTIMIZATION 2: PRE-COMPUTE RF FOR CALIBRATION ===
  calib_qc_preds <- predict(qc_mdl, calib_x, cens_rt)$predictions[,1]

  # === VECTORIZED THRESHOLD EVALUATION ===
  lv1_mat <- predict(mdl,
                     newdata = calib_x,
                     type = "quantile",
                     p = 1 - v_list)

  lv1_mat <- as.matrix(lv1_mat)

  if (nrow(lv1_mat) != nrow(data_calib) && ncol(lv1_mat) == nrow(data_calib)) {
    lv1_mat <- t(lv1_mat)
  }

  qc_mat <- matrix(calib_qc_preds,
                   nrow = nrow(lv1_mat),
                   ncol = ncol(lv1_mat))

  lv_mat <- pmin(lv1_mat, qc_mat)

  if (use_oracle_sc || is_sc_model(mdl0)) {
    pr_mat <- sc_prob(
      mdl0 = mdl0,
      data = data_calib,
      xnames = xnames,
      t = lv_mat
    )
  } else {
    calib_mat <- as.matrix(data_calib[, xnames, drop = FALSE])
    gpr_mean <- mdl0$predict(calib_mat)
    gpr_sd <- mdl0$predict(calib_mat, se.fit = TRUE)$se
    pr_mat <- pnorm((-lv_mat - gpr_mean) / gpr_sd)
  }
  weight_mat <- 1 / pr_mat

  ind1_mat <- (data_calib$censored_T < lv_mat) & (data_calib$C >= lv_mat)
  ind2_mat <- data_calib$C >= lv_mat

  sum_num <- colSums(weight_mat * ind1_mat)
  sum_den <- colSums(weight_mat * ind2_mat)

  alpha_v_list <- sum_num / sum_den

  bad_alpha <- is.na(alpha_v_list) |
    is.nan(alpha_v_list) |
    !is.finite(alpha_v_list) |
    is.na(sum_den) |
    is.nan(sum_den) |
    !is.finite(sum_den) |
    sum_den == 0

  alpha_v_list[bad_alpha] <- 1

  # monotonize alpha
  alpha_v <- monot(alpha_v_list)

  bad_alpha_v <- is.na(alpha_v) |
    is.nan(alpha_v) |
    !is.finite(alpha_v)

  alpha_v[bad_alpha_v] <- 1

  feasible_idx <- which(alpha_v <= alpha)

  if (length(feasible_idx) == 0) {
    v_hat_l <- NULL
  } else {
    v_hat_l <- min(v_list[feasible_idx])
  }

  # === OPTIMIZATION 4: PRE-COMPUTE RF FOR TEST DATA ===
  if (is.null(v_hat_l)) {
    lower_bnd_l <- rep(0, len_x)
  } else {
    newdata_x <- newdata[, xnames, drop = FALSE]
    newdata_qc_preds <- predict(qc_mdl, newdata_x, cens_rt)$predictions[,1]

    # Pass the pre-computed test predictions
    lower_bnd_l <- as.numeric(lv_qct(mdl, newdata, v_hat_l, alpha, newdata_qc_preds))

    # Safety net: If the function somehow still returns a scalar, expand it
    if (length(lower_bnd_l) == 1) {
      lower_bnd_l <- rep(lower_bnd_l, len_x)
    }
  }
  lower_bnd_g <- rep(0, len_x)

  return(list(lower_bnd_l = lower_bnd_l, lower_bnd_g = lower_bnd_g))
}

############################################
## compute lower prediction bound
############################################
lv_qt <- function(mdl, calib_x, v, alpha, cens_rt){
  if(length(v)==0){
    return(lv2_calib = 0)
  }
  lv2_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = 1-v)
  return(lv2_calib)
}

lv_qct <- function(mdl, calib_x, v, alpha, qc_val){
  if(length(v) == 0){
    return(0)
  }
  lv1_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = 1-v)
  lv_calib = pmin(lv1_calib, qc_val)
  return(lv_calib)
}


############################################
## determine the finite candidate set for v
############################################
inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-5)[1]
  }
}

v_pts_qt = function(mdl,
                    data_fit, data_calib,
                    xnames, alpha, cens_rt){

  library(parallel)
  # Check if the operating system is Windows
  if (.Platform$OS.type == "windows") {
    n_threads <- 1  # Windows doesn't support mclapply, force sequential for local testing
  } else {
    # If on Linux/HPC, detect SLURM cores or use all available hardware cores
    slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    n_threads <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
  }

  # Process all rows in parallel
  res_list <- mclapply(1:nrow(data_calib), function(i_calib) {
    # Initialize a vector for this specific row
    pts_row <- c(0, 0) 

    calib_x = data_calib[i_calib, xnames, drop = FALSE]
    names(calib_x) = xnames

    lv = function(v) lv_qt(mdl, calib_x, v, alpha, cens_rt)
    x0 = 0.02
    x1 = 0.95
    upb = lv(x0)
    lwb = lv(x1)
    lv_inv = inverse(function(v) lv_qt(mdl, calib_x, v, alpha, cens_rt), x0, x1)

    if(data_calib$event[i_calib]==0){
      if(data_calib$C[i_calib] >= upb){
        pts_row[1] <- pts_row[2] <- x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts_row[1] <- pts_row[2] <- x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts_row[1] <- pts_row[2] <- lv_inv(data_calib$C[i_calib])$root
      }
    }

    if(data_calib$event[i_calib]==1){
      if(data_calib$C[i_calib] >= upb){
        pts_row[1] <- x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts_row[1] <- x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts_row[1] <- lv_inv(data_calib$C[i_calib])$root
      }
      if(data_calib$censored_T[i_calib] >= upb){
        pts_row[2] <- x0
      }
      if(data_calib$censored_T[i_calib] <= lwb){
        pts_row[2] <- x1
      }
      if((data_calib$censored_T[i_calib] > lwb) && (data_calib$censored_T[i_calib] < upb)){
        pts_row[2] <- lv_inv(data_calib$censored_T[i_calib])$root
      }
    }
    return(pts_row)
  }, mc.cores = n_threads)

  # Bind the list of rows into an N x 2 matrix
  pts <- do.call(rbind, res_list)
  return(pts)
}


v_pts_qct = function(mdl, qc_mdl,
                     data_fit, data_calib,
                     xnames, alpha, cens_rt){

  library(parallel)
  if (.Platform$OS.type == "windows") {
    n_threads <- 1
  } else {
    slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    n_threads <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
  }

  calib_x_full = data_calib[, xnames, drop = FALSE]
  names(calib_x_full) = xnames
  qc_preds <- predict(qc_mdl, calib_x_full, cens_rt)$predictions[,1]

  # Process all rows in parallel
  res_list <- mclapply(1:nrow(data_calib), function(i_calib) {
    pts_row <- c(0, 0)

    calib_x = data_calib[i_calib, xnames, drop = FALSE]
    names(calib_x) = xnames

    # Call the simplified lv_qct
    lv = function(v) lv_qct(mdl, calib_x, v, alpha, qc_preds[i_calib])

    x0 = 0.02
    x1 = 0.95
    upb = lv(x0)
    lwb = lv(x1)

    # FIX: Just pass 'lv' directly into inverse! 
    # This prevents the root finder from calling the Random Forest again.
    lv_inv = inverse(lv, x0, x1)

    if(data_calib$event[i_calib]==0){
      if(data_calib$C[i_calib] >= upb){
        pts_row[1] <- pts_row[2] <- x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts_row[1] <- pts_row[2] <- x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts_row[1] <- pts_row[2] <- lv_inv(data_calib$C[i_calib])$root
      }
    }

    if(data_calib$event[i_calib]==1){
      if(data_calib$C[i_calib] >= upb){
        pts_row[1] <- x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts_row[1] <- x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts_row[1] <- lv_inv(data_calib$C[i_calib])$root
      }
      if(data_calib$censored_T[i_calib] >= upb){
        pts_row[2] <- x0
      }
      if(data_calib$censored_T[i_calib] <= lwb){
        pts_row[2] <- x1
      }
      if((data_calib$censored_T[i_calib] > lwb) && (data_calib$censored_T[i_calib] < upb)){
        pts_row[2] <- lv_inv(data_calib$censored_T[i_calib])$root
      }
    }
    return(pts_row)
  }, mc.cores = n_threads)

  # Bind the list of rows into an N x 2 matrix
  pts <- do.call(rbind, res_list)
  return(pts)
}
