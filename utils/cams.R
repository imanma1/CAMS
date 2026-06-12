cams <- function(x, Xtrain, C, event, time, alpha = 0.1, p, mdl0) {

  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  } else {
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }

  X <- Xtrain
  xnames <- paste0("X", 1:p)
  data <- as.data.frame(cbind(C, event, time, X))
  colnames(data) <- c("C", "event", "censored_T", xnames)

  ## Split the data into the training set and the calibration set
  n <- dim(data)[1]
  n_train <- n/2
  I_fit <- sample(1:n, n_train, replace = FALSE)
  data_fit <- data[I_fit,]
  data_calib <- data[-I_fit,]

  ## Fit the survival model
  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  mdl <- survreg(fmla, data = data_fit, dist = "weibull")
  
  # The truncation level eta for the IPCW weights
  eta = 1/log(n)
  
  start_time = proc.time()[3]

  lower_bnd0 <- est_alpha_ipcw(mdl,
                               newdata[newdata$X1 == 0, , drop=FALSE], 
                               data_calib[data_calib$X1 == 0, , drop=FALSE],
                               xnames, alpha, len_x, mdl0, eta)
                               
  lower_bnd1 <- est_alpha_ipcw(mdl, 
                               newdata[newdata$X1 == 1, , drop=FALSE], 
                               data_calib[data_calib$X1 == 1, , drop=FALSE],
                               xnames, alpha, len_x, mdl0, eta)
                               
  end_time <- proc.time()[3]
  cat(sprintf("est_alpha_ipcw in %.2f seconds.\n", end_time - start_time))

  idx_test_0 <- newdata$X1 == 0
  idx_test_1 <- newdata$X1 == 1
  
  # Initialize a blank numeric vector and fill it
  lower_bnd_vec <- rep(NA, nrow(newdata))
  lower_bnd_vec[idx_test_0] <- lower_bnd0
  lower_bnd_vec[idx_test_1] <- lower_bnd1

  # Convert into a 1-column data frame
  lower_bnd <- data.frame(cams_bnd = lower_bnd_vec)

  return(lower_bnd)
}

# ==========================================
# est_alpha_ipcw
# ==========================================
est_alpha_ipcw <- function(mdl, newdata, data_calib, xnames, alpha, len_x, mdl0, eta) {
  
  # FIX 4: Use a fixed, coarse grid for Lambda as dictated by the 2026 methodology
  v_list <- seq(0.001, 0.999, by = 0.001)
  
  calib_mat <- as.matrix(data_calib[,names(data_calib) %in% xnames, drop=FALSE])
  gpr_mean <- mdl0$predict(calib_mat)
  gpr_sd <- mdl0$predict(calib_mat, se.fit = TRUE)$se

  calib_x <- data_calib[,names(data_calib) %in% xnames, drop=FALSE]
  n_calib_subgroup <- nrow(data_calib)

  pr_calib <- pnorm((gpr_mean - data_calib$censored_T) / gpr_sd)
  pr_calib <- pmax(pr_calib, eta)
  weight_calib <- 1 / pr_calib
  n_eff <- (n_calib_subgroup^2) / sum(weight_calib^2)

  # Pre-calculate effective sample size for the penalty
  # n_eff <- (n_calib_subgroup^2) / sum(weight_calib^2)
  K_maps <- length(v_list)

  # Pre-calculate the Hajek normalization denominator
  # This sums the weights of all actually observed events in the subgroup
  total_weight <- sum(weight_calib[data_calib$event == 1])
  
  est_alpha <- function(v) {
    lv_calib <- lv(mdl, calib_x, v)
    ind = (data_calib$censored_T < lv_calib) & (data_calib$event == 1)
    
    # FIX 1: Hajek Self-Normalization
    sum_num <- sum(weight_calib[ind])
    
    # Prevent division by zero if total_weight is incredibly small
    if(total_weight == 0) {
      risk_empirical <- 1 
    } else {
      risk_empirical <- sum_num / total_weight 
    }
    
    # FIX 2: Reintroduce a Gentle Penalty Buffer
    # We use a small tuning constant C_0 = 0.05 so it doesn't crush the bounds
    C_0 <- 0 
    penalty <- C_0 * sqrt((log(2) + log(K_maps)) / n_eff)
    risk_penalized <- risk_empirical + penalty
    
    if(is.na(risk_penalized)) return(1) else return(risk_penalized)
  }

  alpha_v_list <- unlist(lapply(v_list, est_alpha))
  alpha_v <- cummax(alpha_v_list) # Monotonize
  
  if(sum(alpha_v <= alpha) == 0){
    v_hat_l = NULL
  } else {
    v_hat_l <- max(v_list[alpha_v <= alpha])
  }
  
  if (is.null(v_hat_l)) {
    lower_bnd_l <- rep(0, len_x)
  } else {
    lower_bnd_l <- as.numeric(lv(mdl, newdata, v_hat_l))
    if (length(lower_bnd_l) == 1) {
      lower_bnd_l <- rep(lower_bnd_l, len_x)
    }
  }
  return(lower_bnd_l)
}

# ==========================================
# Core Helper Functions
# ==========================================
lv <- function(mdl, calib_x, v){
  if(length(v) == 0) {
    return(rep(0, times = nrow(calib_x)))
  }
  lv_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = v)
  return(lv_calib)
}