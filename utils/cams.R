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
  
  # 1. Dynamically calculate v_list based on exact crossing points
  pts <- v_pts_cams(mdl, data_calib, xnames)
  v_list <- sort(unique(c(as.vector(pts), 0.02, 0.95))) # Flatten, sort, and cap at bounds
  
  calib_mat <- as.matrix(data_calib[,names(data_calib) %in% xnames, drop=FALSE])
  gpr_mean <- mdl0$predict(calib_mat)
  gpr_sd <- mdl0$predict(calib_mat, se.fit = TRUE)$se

  calib_x <- data_calib[,names(data_calib) %in% xnames, drop=FALSE]
  n_calib_subgroup <- nrow(data_calib) 

  est_alpha <- function(v) {
    lv_calib <- lv(mdl, calib_x, v)
    
    pr_calib <- pnorm((-lv_calib - gpr_mean) / gpr_sd)
    pr_calib <- pmax(pr_calib, eta) 
    weight_calib <- 1 / pr_calib
    
    ind = (data_calib$censored_T < lv_calib) & (data_calib$event == 1)
    
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

  est_start_time <- proc.time()[3]
  alpha_v_list <- unlist(mclapply(v_list, est_alpha, mc.cores = n_threads))
  cat(sprintf("est_alpha in %.2f seconds\n", proc.time()[3] - est_start_time))
  
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

# ==========================================
# Core Helper Functions
# ==========================================
lv <- function(mdl, calib_x, v){
  if(length(v) == 0) {
    return(rep(0, times = nrow(calib_x)))
  }
  lv_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = 1 - v)
  return(lv_calib)
}

inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-5)[1]
  }
}

v_pts_cams = function(mdl, data_calib, xnames){
  library(parallel)
  if (.Platform$OS.type == "windows") {
    n_threads <- 1
  } else {
    slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    n_threads <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
  }
  
  res_list <- mclapply(1:nrow(data_calib), function(i_calib) {
    pts_row <- c(0, 0) 
    
    calib_x = data_calib[i_calib, names(data_calib) %in% xnames, drop=FALSE]
    names(calib_x) = xnames
    
    lv_local = function(v) lv(mdl, calib_x, v)
    x0 = 0.02
    x1 = 0.95
    upb = lv_local(x0)
    lwb = lv_local(x1)
    
    lv_inv = inverse(lv_local, x0, x1)
    
    C_val = data_calib$C[i_calib]
    T_val = data_calib$censored_T[i_calib]
    evt = data_calib$event[i_calib]
    
    if(evt == 0){
      if(C_val >= upb) pts_row[1] <- pts_row[2] <- x0
      if(C_val <= lwb) pts_row[1] <- pts_row[2] <- x1
      if(C_val > lwb && C_val < upb) pts_row[1] <- pts_row[2] <- lv_inv(C_val)$root
    } else {
      if(C_val >= upb) pts_row[1] <- x0
      if(C_val <= lwb) pts_row[1] <- x1
      if(C_val > lwb && C_val < upb) pts_row[1] <- lv_inv(C_val)$root
      
      if(T_val >= upb) pts_row[2] <- x0
      if(T_val <= lwb) pts_row[2] <- x1
      if(T_val > lwb && T_val < upb) pts_row[2] <- lv_inv(T_val)$root
    }
    return(pts_row)
  }, mc.cores = n_threads)
  
  pts <- do.call(rbind, res_list)
  return(pts)
}