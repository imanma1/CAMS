cams <- function(x,
  Xtrain, C, event, time,
  alpha = 0.1, p,
  mdl0
){

  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
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
  fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
                           paste(xnames, collapse= "+")))
  mdl <- survreg(fmla, data = data_fit, dist = "weibull")
  
  # The truncation level eta for the IPCW weights
  eta = 1/log(n)
  
  lower_bnd0 <- est_alpha_ipcw(mdl,
                               newdata[newdata$X1 == 0, , drop=FALSE],
                               data_calib[data_calib$X1 == 0, , drop=FALSE],
                               xnames, alpha, len_x, mdl0, eta)
  lower_bnd1 <- est_alpha_ipcw(mdl,
                               newdata[newdata$X1 == 1, , drop=FALSE],
                               data_calib[data_calib$X1 == 1, , drop=FALSE],
                               xnames, alpha, len_x, mdl0, eta)

  idx_test_0 <- newdata$X1 == 0
  idx_test_1 <- newdata$X1 == 1
  lower_bnd_vec <- rep(NA, nrow(newdata))
  lower_bnd_vec[idx_test_0] <- lower_bnd0
  lower_bnd_vec[idx_test_1] <- lower_bnd1
  lower_bnd <- data.frame(cams_bnd = lower_bnd_vec)

  return(lower_bnd)
}

est_alpha_ipcw <- function(mdl, newdata, data_calib,
                  xnames, alpha, len_x, mdl0, eta) {
  v_list <- seq(0, 1, by = 0.01)

  calib_mat <- as.matrix(data_calib[,names(data_calib) %in% xnames, drop=FALSE])
  gpr_mean <- mdl0$predict(calib_mat)
  gpr_sd <- mdl0$predict(calib_mat, se.fit = TRUE)$se

  calib_x <- data_calib[,names(data_calib) %in% xnames, drop=FALSE]
  n_calib_subgroup <- nrow(data_calib) # Subgroup size n_R

  est_alpha <- function(v) {
    lv_calib <- lv(mdl, calib_x, v)
    
    # Calculate censoring probabilities and apply the eta truncation
    pr_calib <- pnorm((-lv_calib - gpr_mean) / gpr_sd)
    pr_calib <- pmax(pr_calib, eta) 
    weight_calib <- 1 / pr_calib
    
    # The IPCW indicator MUST include the event status (Delta = 1)
    ind = (data_calib$censored_T < lv_calib) & (data_calib$event == 1)
    
    # Sum the weights and divide by the subgroup sample size to get the risk
    sum_num <- sum(weight_calib[ind])
    risk <- sum_num / n_calib_subgroup
    
    if(is.na(risk)) return(1) else return(risk)
  }

  library(parallel)
  if (.Platform$OS.type == "windows") {
    n_threads <- 1
  } else {
    slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    n_threads <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
  }
  alpha_v_list <- unlist(mclapply(v_list, est_alpha, mc.cores = n_threads))
  
  # monotonize alpha
  alpha_v <- monot(alpha_v_list)
  if(sum(alpha_v <= alpha) == 0){
    v_hat_l = NULL
  }else{
    v_hat_l <- min(v_list[alpha_v <= alpha])
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

lv <- function(mdl, calib_x, v){
  # Simplified to just return the quantile prediction
  if(length(v) == 0) {
    return(rep(0, times = nrow(calib_x)))
  }
  lv_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = 1 - v)
  return(lv_calib)
}