simu <- function(seed, setting, n, p,
                 n_train, n_calib, n_test,
                 beta, xmin, xmax,
                 exp_rate, alpha){
  mod <- "cox"

  ## Initialization
  set.seed(seed)
  xnames <- paste0("X", 1:p)

  ## Generate data according to the setting
  data_obj <- model_generating_fun(n_train, n_calib, n_test,
                                   setting, beta, xnames,
                                   xmin, xmax, exp_rate)
  p <- p + 1 # account for the added bernoulli variable
  xnames <- paste0("X", 1:p)
  xnames_sub <- paste0("X", 2:p) # Feature names without X1 for subgroup training

  data_fit <- data_obj$data_fit
  data_calib <- data_obj$data_calib
  data_test <- data_obj$data_test
  data <- data_obj$data
  T_test <- data_obj$T_test
  
  alpha_list <- alpha
  
  ########################################
  ## Core Pipeline Helper Function
  ########################################
  # This function trains all 6 methods and returns their lower bounds and times
  run_pipeline <- function(sub_fit, sub_calib, sub_test, sub_data, xnames_to_use, alpha_list, seed, mod) {
    # Utility function to safely extract quantiles from a coxph object
    extract_quant <- function(mdl, x, alpha){
      res <- summary(survfit(mdl, newdata = x))
      time_point <- res$time
      survcdf <- 1 - res$surv
      if(sum(survcdf >= alpha)==0){
        quant = max(time_point)
      }else{
        quant <- time_point[min(which(survcdf >= alpha))]
      }
      return(quant)
    }
    
    n_train_sub <- nrow(sub_fit)
    p_sub <- length(xnames_to_use)
    
    x <- sub_test[, xnames_to_use, drop=FALSE]
    Xtrain <- sub_data[, xnames_to_use, drop=FALSE]
    C <- sub_data$C
    event <- sub_data$event
    time <- sub_data$censored_T
    
    fit <- sub_fit
    fit$C <- -fit$C
    
    # 1. Estimate mdl0
    start_time <- proc.time()[3]
    mdl0 <- GauPro(X = as.matrix(fit[, xnames_to_use, drop=FALSE]), Z = fit$C, D = p_sub, type = "Gauss")
    time_mdl0 <- proc.time()[3] - start_time
    
    # 2. cfsurv_c (qt and qct)
    start_time <- proc.time()[3]
    lb_res <- cfsurv_c(x=x, Xtrain=Xtrain, C=C, event=event, time=time, alpha=alpha_list, mdl0=mdl0)
    time_q <- proc.time()[3] - start_time
    
    res1 <- lb_res$lower_bnd_qtl
    res2 <- lb_res$lower_bnd_qctl
    time_qt <- lb_res$time_qt
    time_qct <- lb_res$time_qct
    time_q_base <- time_q - (time_qt + time_qct)
    time_qt <- time_qt + time_q_base + time_mdl0
    time_qct <- time_qct + time_q_base + time_mdl0
    
    times <- c(time_qt, time_qct)
    output <- data.frame(qtl = res1, qctl = res2)
    
    # 3. cfsurv (qc0)
    start_time <- proc.time()[3]
    res0 <- cfsurv(x, c_list=NULL, pr_list=NULL, pr_new_list=NULL,
                   Xtrain, C, event, time, alpha=alpha_list, type="quantile",
                   model=mod, dist="weibull", I_fit=NULL, ftol=.1, tol=.1,
                   n.tree=100, mdl0=mdl0)
    time_qc0 <- proc.time()[3] - start_time + time_mdl0
    times <- c(times, time_qc0)
    output$qc0 <- res0$res
    
    # 4. vanilla CQR
    start_time <- proc.time()[3]
    res <- lapply(alpha_list, cqr,
                  x = x,
                  Xtrain = Xtrain,
                  Ytrain = sub_data$censored_T,
                  I_fit = 1:n_train_sub,
                  seed = seed + 7)
    res <- do.call(rbind, lapply(res, as.data.frame))
    output$cqr.bnd <- res[, 1]
    times <- c(times, proc.time()[3] - start_time)
    
    # 5. Cox Model
    start_time <- proc.time()[3]
    fmla <- as.formula(paste("Surv(censored_T, event) ~", paste(xnames_to_use, collapse="+")))
    mdl <- coxph(fmla, data = sub_data)
    cox_res <- c()
    for (i in 1:nrow(sub_test)) {
       cox_res <- c(cox_res, extract_quant(mdl, sub_test[i, , drop=FALSE], alpha_list))
    }
    output$cox.bnd <- cox_res
    times <- c(times, proc.time()[3] - start_time)
    
    # 6. Random Forest
    start_time <- proc.time()[3]
    ntree <- 1000
    nodesize <- 80
    fmla_rf <- as.formula(paste("censored_T ~", paste(xnames_to_use, collapse="+")))
    mdl <- crf.km(fmla_rf, ntree = ntree, nodesize = nodesize,
                  data_train = sub_data[, c(xnames_to_use, "censored_T", "event"), drop=FALSE], 
                  data_test = sub_test[, xnames_to_use, drop=FALSE], 
                  yname = 'censored_T', iname = 'event', tau = alpha_list, method = "grf")
    output$rf.bnd <- mdl$predicted
    times <- c(times, proc.time()[3] - start_time)
    
    return(list(output = output, times = times))
  }
  
  ########################################
  ## Track Subgroup Indices
  ########################################
  idx_test_0 <- data_test$X1 == 0
  idx_test_1 <- data_test$X1 == 1

  ########################################
  ## APPROACH 1: Joint Modeling
  ########################################
  cat("========== Executing Approach 1: Joint Modeling ==========\n")
  res_joint <- run_pipeline(data_fit, data_calib, data_test, data, 
                            xnames, alpha_list, seed, mod)

  ########################################
  ## APPROACH 2: Subgroup Modeling
  ########################################
  cat("========== Executing Approach 2: Subgroup X1 = 0 ==========\n")
  res_0 <- run_pipeline(data_fit[data_fit$X1 == 0, , drop=FALSE],
                        data_calib[data_calib$X1 == 0, , drop=FALSE],
                        data_test[data_test$X1 == 0, , drop=FALSE],
                        data[data$X1 == 0, , drop=FALSE],
                        xnames_sub, alpha_list, seed, mod)
                                 
  cat("========== Executing Approach 2: Subgroup X1 = 1 ==========\n")
  res_1 <- run_pipeline(data_fit[data_fit$X1 == 1, , drop=FALSE],
                        data_calib[data_calib$X1 == 1, , drop=FALSE],
                        data_test[data_test$X1 == 1, , drop=FALSE],
                        data[data$X1 == 1, , drop=FALSE],
                        xnames_sub, alpha_list, seed, mod)
  
  # Merge subgroup outputs to map exactly to the data_test row order
  output_subgroup <- data.frame(matrix(ncol = ncol(res_0$output), nrow = nrow(data_test)))
  colnames(output_subgroup) <- colnames(res_0$output)
  output_subgroup[idx_test_0, ] <- res_0$output
  output_subgroup[idx_test_1, ] <- res_1$output
  times_subgroup <- res_0$times + res_1$times

  ########################################
  ## Format Output Helper
  ########################################
  compute_metrics <- function(output_df, times_vec, suffix_label) {
    # Coverage logic
    cov_marg <- apply(output_df, 2, function(x) sum(T_test >= x) / length(x))
    cov_grp0 <- apply(output_df, 2, function(x) sum(T_test[idx_test_0] >= x[idx_test_0]) / sum(idx_test_0))
    cov_grp1 <- apply(output_df, 2, function(x) sum(T_test[idx_test_1] >= x[idx_test_1]) / sum(idx_test_1))
    
    # Mean logic
    simulen      <- apply(output_df, 2, mean)
    simulen_grp0 <- apply(output_df, 2, function(x) mean(x[idx_test_0]))
    simulen_grp1 <- apply(output_df, 2, function(x) mean(x[idx_test_1]))
    
    method_names <- paste(c("DFT-adaptive-T", "DFT-adaptive-CT", "DFT-fixed", 
                            "Vanilla CQR", "Cox", "Random Forest"), suffix_label)
    
    data.frame(
      "method"                       = method_names,
      "setting"                      = setting,
      "group coverage for x_1 = 0"   = cov_grp0,
      "group coverage for x_1 = 1"   = cov_grp1,
      "Marginal coverage"            = cov_marg,
      "lower bound mean for x_1 = 0" = simulen_grp0,
      "lower bound mean for x_1 = 1" = simulen_grp1,
      "lower bound values mean"      = simulen,
      "computation time"             = times_vec,
      check.names = FALSE
    )
  }
  
  ########################################
  ## Compute & Bind Final Results
  ########################################
  df_joint <- compute_metrics(res_joint$output, res_joint$times, "(Joint)")
  df_subgroup <- compute_metrics(output_subgroup, times_subgroup, "(Subgroup)")
  
  simu_out <- rbind(df_joint, df_subgroup)
  rownames(simu_out) <- NULL
  
  return(simu_out)
}