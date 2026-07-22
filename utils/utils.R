#' Fitting the censoring probability P(C>=c|X)
#'
#' @export

censoring_prob <- function(fit, calib, test=NULL,
                           method="gpr",
                           xnames, c,
                           ftol=.1, tol=.1, n.tree = 40){

  p <- length(xnames)
  if(method == "np"){
    ## Fitting P(-C<=-c_0|X) (since P(C>=c_0|X)=P(-C<=-c_0|X))
    fit$C <- -fit$C
    fmla <- with(fit,as.formula(paste("C ~ ", paste(xnames, collapse= "+"))))
    if(length(xnames)==1){
      capture.output(bw <- npcdistbw(fmla),file =NULL)
    }else{
      capture.output(bw <- npcdistbw(fmla,ftol=ftol,tol=tol),file =NULL)
    }

    ## Computing censoring scores for the calibration data
    newdata_calib <- calib
    newdata_calib$C <- -c
    pr_calib<- npcdist(bws=bw,newdata = newdata_calib)$condist

    ## Computing the censoring scores for the test data
    if(!is.null(test)){
      newdata <- cbind(test,C=-c)
      newdata <- data.frame(newdata)
      colnames(newdata) <- c(xnames,"C")
      pr_new <- npcdist(bws=bw,newdata=newdata)$condist
    }else{pr_new=NULL}
  }

  if(method == "distBoost"){
    ## Fitting P(-C<=-c_0|X) (since P(C>=c_0|X)=P(-C<=-c_0|X))
    fit$C <- -fit$C
    fmla <- with(fit,as.formula(paste("C ~ ", paste(xnames, collapse= "+"))))
    gbm_mdl <- gbm(fmla,data=fit,distribution="gaussian", n.tree = n.tree)
    median_fit<- predict(object=gbm_mdl,newdata = fit)
    res_fit <- fit$C-median_fit
    resamp_fit <- median_fit + res_fit[sample.int(dim(fit)[1])]
    if(p==1){
      xdf <- data.frame(X1=fit[,colnames(fit)%in%xnames],
                    X2=rep(1,dim(fit)[1]))
     }else{
      xdf <- fit[,colnames(fit)%in%xnames]
     }
    mdlrb <- modtrast(xdf,fit$C,resamp_fit,min.node=200)
    
    ## Computing the censoring scores for the fitting data
    pr_fit <- rep(NA, dim(fit)[1])
    for(i in 1:length(pr_fit)){
        pr_fit[i] <- distBoost_cdf(mdlrb,xdf[i,],median_fit[i],-c[i],res_fit)
    }

    ## Computing the censoring scores for the calibration data
    pr_calib <- rep(NA, dim(calib)[1])
    median_calib<- predict(object=gbm_mdl,newdata = calib)
    if(p==1){
      xdf <- data.frame(X1=calib[,colnames(calib)%in%xnames],
                    X2=rep(1,dim(calib)[1]))
     }else{
      xdf <- calib[,colnames(calib)%in%xnames]
     }

    for(i in 1:length(pr_calib)){
      pr_calib[i] <- distBoost_cdf(mdlrb,xdf[i,],median_calib[i],
                               -c[i],res_fit)
    }

    ## Computing the censoring scores for the test data
    if(!is.null(test)){
      newdata <- data.frame(test)
      colnames(newdata) <- xnames
      median_test<- predict(object=gbm_mdl,newdata = newdata)
      if(p==1){
        xdf <- data.frame(X1=test,
                    X2=rep(1,length(test)))
        n_new <- length(test)
      }else{
        xdf <- test
        n_new <- dim(test)[1]
      }
      pr_new <- rep(NA, n_new)
      for(i in 1:n_new){
          pr_new[i] <- distBoost_cdf(mdlrb,xdf[i,],median_test[i],-c,res_fit)
      }
    }else{pr_new=NULL}
  }


  if(method == "gpr"){
    ## Fitting P(-C<=-c_0|X) (since P(C>=c_0|X)=P(-C<=-c_0|X))
    fit$C <- -fit$C
    gpr_mdl <- GauPro(X = as.matrix(fit[,names(fit) %in% xnames]),
                    Z = fit$C, D = p,
                    type = "Gauss")
    
    ## Computing the censoring scores for the fitting data
    mean_fit <- gpr_mdl$predict(as.matrix(fit[,names(fit) %in% xnames]))
    sd_fit <- gpr_mdl$predict(as.matrix(fit[,names(fit) %in% xnames]),
                              se.fit = TRUE)$se

    pr_fit <- pnorm((-c - mean_fit) / sd_fit)

    ## Computing the censoring scores for the calibration data
    mean_calib <- gpr_mdl$predict(as.matrix(calib[,names(calib) %in% xnames]))
    sd_calib <- gpr_mdl$predict(as.matrix(calib[,names(calib) %in% xnames]),
                              se.fit = TRUE)$se

    pr_calib <- pnorm((-c - mean_calib) / sd_calib)

    ## Computing the censoring scores for the test data
    if(!is.null(test)){
      newdata <- data.frame(test)
      colnames(newdata) <- xnames
      mean_new <- gpr_mdl$predict(as.matrix(newdata[,names(newdata) %in% xnames]))
      sd_new <- gpr_mdl$predict(as.matrix(newdata[,names(newdata) %in% xnames]),
                              se.fit = TRUE)$se

      pr_new <- pnorm((-c - mean_new) / sd_new)

    }else{pr_new=NULL}

  }
  return(list(pr_fit = pr_fit, pr_calib = pr_calib, pr_new = pr_new))
}

cens_prob <- function(mdl, calib, test=NULL,
                      method="gpr",
                      xnames, c,
                      ftol=.1, tol=.1, n.tree = 40){
  
  p <- length(xnames)
  gpr_mdl = mdl
  ## Computing the censoring scores for the calibration data
  mean_calib <- gpr_mdl$predict(as.matrix(calib[,names(calib) %in% xnames]))
  sd_calib <- gpr_mdl$predict(as.matrix(calib[,names(calib) %in% xnames]),
                              se.fit = TRUE)$se
  
  pr_calib <- pnorm((-c - mean_calib) / sd_calib)
  
  ## Computing the censoring scores for the test data
  if(!is.null(test)){
    newdata <- data.frame(test)
    colnames(newdata) <- xnames
    mean_new <- gpr_mdl$predict(newdata[,names(newdata) %in% xnames])
    sd_new <- gpr_mdl$predict(newdata[,names(newdata) %in% xnames],
                              se.fit = TRUE)$se
    
    pr_new <- pnorm((-c - mean_new) / sd_new)
    
  }else{pr_new=NULL}
  return(list(pr_calib = pr_calib, pr_new = pr_new))
}


#' Computing the calibration term with covaraite shift
#'
#' construct the one-sided confidence interval for a unit's survival time T
#'
#' @param x a vector of the covariate of the test data.
#' @param r the censoring time of the test data.
#' @param alpha a number betweeo 0 and 1, specifying the miscaverage rate.
#' @param data a data frame used for calibration, containing four columns: (X,R,event,censored_T). 
#' @param mdl The fitted model to estimate the conditional quantile (default is NULL).
#' @param quant_lo the fitted conditional quantile for the calibration data (default is NULL).
#' @param new_quant_lo the fitted conditional quantile for the test data (default is NULL).
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @family confint
#'
#' @export


get_calibration <- function(score,weight_calib,weight_new,alpha){
  ## Check input format
  if(length(score)!=length(weight_calib)) stop("The length of score is not compatible with the length of weight!")

  if(!is.numeric(alpha)) stop("alpha should be a real number between 0 and 1!")
  if(alpha>1 | alpha<0) stop("alpha should be a real number between 0 and 1!")

  ## Computing the calibration term
  weight <- c(weight_calib,weight_new)
  weight <- weight/sum(weight)
  score_vec <- c(score,Inf)
  sort_score <- sort(score_vec)
  order_score <- order(score_vec)
  sort_w <- weight[order_score]
  idxw <- min(which(cumsum(sort_w)>=1-alpha))
  calib_term <- sort_score[idxw]

  return(calib_term)
}

 

get_survival_fun <- function(x,t,bw,xnames){
  input_data <- data.frame(x)
  colnames(input_data) <- xnames
  input_data <- cbind(input_data,censored_T=t)
  val<- npcdist(bws=bw,newdata=input_data)$condist
  return(val)
}



monot = function(a){
  m_a = NULL
  for(i in 1:length(a)){
    m_a = c(m_a, max(a[i:length(a)]))
  }
  return(m_a)
}


inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-5)[1]
  }
}

make_oracle_sc_model <- function(setting) {
  structure(
    list(
      setting = setting,
      sigma_c = 0.5
    ),
    class = "oracle_sc"
  )
}

make_power_oracle_sc_model <- function(setting, gamma) {

  if (!is.finite(gamma) || gamma <= 0) {
    stop("gamma must be positive.")
  }

  structure(
    list(
      setting = setting,
      sigma_c = 0.5,
      gamma = gamma
    ),
    class = "power_oracle_sc"
  )
}


mu_c_oracle <- function(X, setting) {
  X <- as.data.frame(X)

  if (setting == "homo_cens") {
    mu_c <- rep(3.0, nrow(X))

  } else if (setting == "cov_cens") {
    mu_c <- 2.5 + 1.2 * X$X2

  } else if (setting == "prot_cens") {
    mu_c <- 3.0 - 1.5 * X$X1

  } else if (setting == "heavy_prot_cens") {
    mu_c <- 3.5 - 2.8 * X$X1

  } else if (setting == "heavy_inter_cens") {
    mu_c <- 3.5 - 2.8 * ((X$X1 == 1) & (X$X2 > 0))

  } else if (setting == "surv_misspec") {
    mu_c <- 2.5 + 0.5 * X$X2

  } else if (setting == "cens_misspec") {
    mu_c <- 2.5 +
      (X$X2^2) * X$X1 +
      cos(X$X3) +
      0.2 * rowSums(X[, paste0("X", 4:10), drop = FALSE]^2)

  } else if (setting == "simul_misspec") {
    mu_c <- 2.5 +
      (X$X3^2) -
      1.5 * X$X1 +
      0.1 * rowSums(abs(X[, paste0("X", 4:10), drop = FALSE]))

  } else if (setting == "complex_surv") {
    mu_c <- rep(3.0, nrow(X))

  } else if (setting == "var_shift_heavy_cens") {
    mu_c <- 3.0 + 0.5 * X$X2 - 2.8 * X$X1

  } else if (setting == "starve_hetero") {
    mu_c <- 3.0 - 1.5 * X$X1

  } else if (setting == "starve_hetero_high_dim") {
    mu_c <- 3.0 - 1.5 * X$X1

  } else if (setting == "high_survival_heavy_cens") {
    mu_c <- 3.40 -
      0.50 * X$X2 +
      0.20 * X$X3 -
      0.80 * X$X1

  } else if (setting == "anti_aligned_cens") {
    mu_c <- 2.80 -
      0.70 * X$X2 +
      0.50 * X$X3 -
      0.80 * X$X1

  } else if (setting == "weibull_aft_anti_cens") {
    mu_c <- 3.30 -
      0.60 * X$X2 +
      0.40 * X$X3 -
      0.80 * X$X1

  } else if (setting == "moderate_inter_cens") {
    mu_c <- 3.00 -
      0.40 * ((X$X2 > 0) & (X$X3 < 0)) -
      0.60 * X$X1 +
      0.20 * X$X4
  } else if (setting == "cams_pool_calib_hd") {

    p <- 75

    beta <- 0.05 * rep(c(1, -1), length.out = p - 1)

    score <- as.numeric(
      as.matrix(
        X[, paste0("X", 2:p), drop = FALSE]
      ) %*% beta
    )

    mu_c <- 3.15 -
      0.15 * X$X1 -
      0.10 * score
  } else if (setting == "cams_pool_calib_hd_mild") {

    p <- 75
    beta <- 0.05 * rep(c(1, -1), length.out = p - 1)

    score <- as.numeric(
      as.matrix(
        X[, paste0("X", 2:p), drop = FALSE]
      ) %*% beta
    )

    mu_c <- 3.20 -
      0.10 * X$X1 -
      0.08 * score
  } else if (setting == "cams_pool_calib_hd_strong") {

    p <- 75
    beta <- 0.05 * rep(c(1, -1), length.out = p - 1)

    score <- as.numeric(
      as.matrix(
        X[, paste0("X", 2:p), drop = FALSE]
      ) %*% beta
    )

    mu_c <- 3.20 -
      0.10 * X$X1 -
      0.10 * score
  } else if (
    setting == "rare_intersection_shared_weibull"
  ) {

    group_1_positive <- as.numeric(
      X$X1 == 1 &
        X$X2 > 0
    )

    group_1_nonpositive <- as.numeric(
      X$X1 == 1 &
        X$X2 <= 0
    )

    group_0_positive <- as.numeric(
      X$X1 == 0 &
        X$X2 > 0
    )

    mu_c <- 3.20 +
      0.20 * X$X3 -
      0.60 * group_1_positive -
      0.30 * group_1_nonpositive -
      0.20 * group_0_positive

  } else if (
    setting == "basis_intersection_shared_weibull"
  ) {

    group_1_positive <- as.numeric(
      X$X1 == 1 &
        X$X2 > 0
    )

    group_0_positive <- as.numeric(
      X$X1 == 0 &
        X$X2 > 0
    )

    mu_c <- 3.00 -
      0.65 * group_1_positive -
      0.35 * group_0_positive +
      0.25 * X$X3
  } else if (
    setting %in% c(
      "intersection_scale_shift_hd",
      "intersection_early_event_mixture_hd"
    )
  ) {

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
          "Oracle censoring model is missing: %s",
          paste(
            missing_names,
            collapse = ", "
          )
        )
      )
    }

    mu_c <- 3.25 -
      0.15 * X$X1 -
      0.10 * as.numeric(X$X2 > 0) +
      0.15 * X$X3
  } else {
    stop(sprintf("Unknown setting for oracle S_C: %s", setting))
  }

  as.numeric(mu_c)
}

two_component_lognormal_survival <- function(
    t,
    mixture_probability,
    first_mean,
    first_sd,
    second_mean,
    second_sd
) {

  n <- length(mixture_probability)

  if (
    length(first_mean) != n ||
      length(second_mean) != n
  ) {
    stop(
      "Mixture probabilities and component means must have equal lengths."
    )
  }

  if (is.matrix(t)) {

    if (nrow(t) != n) {
      stop(
        paste0(
          "For matrix t, nrow(t) must equal ",
          "the number of observations."
        )
      )
    }

    k <- ncol(t)

    probability_matrix <- matrix(
      mixture_probability,
      nrow = n,
      ncol = k
    )

    first_mean_matrix <- matrix(
      first_mean,
      nrow = n,
      ncol = k
    )

    second_mean_matrix <- matrix(
      second_mean,
      nrow = n,
      ncol = k
    )

    log_t <- log(
      pmax(
        t,
        .Machine$double.xmin
      )
    )

    first_survival <- 1 - pnorm(
      (
        log_t -
          first_mean_matrix
      ) / first_sd
    )

    second_survival <- 1 - pnorm(
      (
        log_t -
          second_mean_matrix
      ) / second_sd
    )

    survival_probability <- (
      probability_matrix * first_survival +
        (
          1 - probability_matrix
        ) * second_survival
    )

    survival_probability[t <= 0] <- 1

  } else {

    if (length(t) == 1L) {
      t <- rep(t, n)
    }

    if (length(t) != n) {
      stop(
        paste0(
          "Length of t must equal the ",
          "number of observations."
        )
      )
    }

    log_t <- log(
      pmax(
        t,
        .Machine$double.xmin
      )
    )

    first_survival <- 1 - pnorm(
      (
        log_t -
          first_mean
      ) / first_sd
    )

    second_survival <- 1 - pnorm(
      (
        log_t -
          second_mean
      ) / second_sd
    )

    survival_probability <- (
      mixture_probability * first_survival +
        (
          1 - mixture_probability
        ) * second_survival
    )

    survival_probability[t <= 0] <- 1
  }

  survival_probability[
    !is.finite(survival_probability)
  ] <- 0

  pmin(
    pmax(
      survival_probability,
      0
    ),
    1
  )
}

oracle_sc_prob <- function(oracle_mdl, data, t) {

  X <- as.data.frame(data)
  setting <- oracle_mdl$setting

  if (
    setting == "mixture_intersection_shared_weibull"
  ) {

    required_names <- paste0(
      "X",
      1:20
    )

    missing_names <- setdiff(
      required_names,
      colnames(X)
    )

    if (length(missing_names) > 0L) {
      stop(
        sprintf(
          paste0(
            "Oracle mixture censoring model ",
            "is missing columns: %s"
          ),
          paste(
            missing_names,
            collapse = ", "
          )
        )
      )
    }

    beta_dense <- 0.05 * (-1)^(5:20)

    dense_score <- as.numeric(
      as.matrix(
        X[
          ,
          paste0("X", 5:20),
          drop = FALSE
        ]
      ) %*% beta_dense
    )

    mu_t <- 2.40 +
      0.50 * X$X1 +
      0.70 * X$X2 -
      0.50 * X$X3 +
      0.30 * X$X4 +
      dense_score

    group_1_positive <- as.numeric(
      X$X1 == 1 &
        X$X2 > 0
    )

    group_0_positive <- as.numeric(
      X$X1 == 0 &
        X$X2 > 0
    )

    prob_early <- 0.08 +
      0.25 * group_1_positive +
      0.12 * group_0_positive

    prob_early <- pmin(
      pmax(prob_early, 0),
      1
    )

    return(
      two_component_lognormal_survival(
        t = t,
        mixture_probability = prob_early,
        first_mean = mu_t - 0.35,
        first_sd = 0.30,
        second_mean = mu_t + 1.00,
        second_sd = 0.40
      )
    )
  }

  mixture_settings <- c(
    "cams_vs_vanilla_lower_tail_hd_mild",
    "cams_vs_vanilla_lower_tail_hd_main",
    "cams_vs_vanilla_lower_tail_hd_strong",
    "local_null_constant_scale",
    "local_within_group_scale_x2",
    "local_interaction_scale_x2_x3",
    "local_unsupported_rare_pocket",
    "local_oblique_scale_x2_x3"
  )

  # ==========================================================
  # Special oracle for the two-component censoring mixtures
  # ==========================================================
  if (setting %in% mixture_settings) {

    p <- 75

    required_names <- paste0("X", 1:p)

    missing_names <- setdiff(required_names, colnames(X))

    if (length(missing_names) > 0) {
      stop(
        sprintf(
          "Oracle mixture censoring model is missing columns: %s",
          paste(missing_names, collapse = ", ")
        )
      )
    }

    beta_dense <- 0.03 * rep(
      c(1, -1),
      length.out = p - 4
    )

    dense_score <- as.numeric(
      as.matrix(
        X[, paste0("X", 5:p), drop = FALSE]
      ) %*% beta_dense
    )

    mu_t <- 2.8 +
      0.40 * X$X1 +
      0.60 * X$X2 -
      0.50 * X$X3 +
      0.30 * X$X4 +
      dense_score

  if (setting == "cams_vs_vanilla_lower_tail_hd_mild") {

    prob_early <- 0.10 + 0.05 * X$X1
    early_offset <- 1.00

  } else if (
    setting %in% c(
      "cams_vs_vanilla_lower_tail_hd_main",
      "local_null_constant_scale",
      "local_within_group_scale_x2",
      "local_interaction_scale_x2_x3",
      "local_unsupported_rare_pocket",
      "local_oblique_scale_x2_x3"
    )
  ) {

    # All Local-CAMS settings use the main censoring mechanism
    prob_early <- 0.15 + 0.05 * X$X1
    early_offset <- 1.30

  } else if (setting == "cams_vs_vanilla_lower_tail_hd_strong") {

    prob_early <- 0.20 + 0.05 * X$X1
    early_offset <- 1.50

  } else {

    stop(
      sprintf(
        "Unknown mixture-censoring setting: %s",
        setting
      )
    )
  }

    mu_early <- mu_t - early_offset
    mu_late <- mu_t + 1.20

    sigma_early <- 0.25
    sigma_late <- 0.40

    if (is.matrix(t)) {

      if (nrow(t) != nrow(X)) {
        stop(
          sprintf(
            "For matrix t, nrow(t) must equal nrow(data): %d versus %d.",
            nrow(t),
            nrow(X)
          )
        )
      }

      n <- nrow(X)
      k <- ncol(t)

      prob_early_mat <- matrix(
        prob_early,
        nrow = n,
        ncol = k
      )

      mu_early_mat <- matrix(
        mu_early,
        nrow = n,
        ncol = k
      )

      mu_late_mat <- matrix(
        mu_late,
        nrow = n,
        ncol = k
      )

      t_safe <- pmax(t, .Machine$double.xmin)
      log_t <- log(t_safe)

      surv_early <- 1 - pnorm(
        (log_t - mu_early_mat) / sigma_early
      )

      surv_late <- 1 - pnorm(
        (log_t - mu_late_mat) / sigma_late
      )

      pr <- prob_early_mat * surv_early +
        (1 - prob_early_mat) * surv_late

      pr[t <= 0] <- 1

    } else {

      if (length(t) == 1) {
        t <- rep(t, nrow(X))
      }

      if (length(t) != nrow(X)) {
        stop(
          sprintf(
            "Length of t must equal nrow(data): %d versus %d.",
            length(t),
            nrow(X)
          )
        )
      }

      t_safe <- pmax(t, .Machine$double.xmin)
      log_t <- log(t_safe)

      surv_early <- 1 - pnorm(
        (log_t - mu_early) / sigma_early
      )

      surv_late <- 1 - pnorm(
        (log_t - mu_late) / sigma_late
      )

      pr <- prob_early * surv_early +
        (1 - prob_early) * surv_late

      pr[t <= 0] <- 1
    }

    pr[!is.finite(pr)] <- 0

    return(
      pmin(
        pmax(pr, 0),
        1
      )
    )
  }

  # ==========================================================
  # Existing single-lognormal oracle settings
  # ==========================================================
  mu_c <- mu_c_oracle(X, setting)
  sigma_c <- oracle_mdl$sigma_c

  if (is.matrix(t)) {

    if (nrow(t) != nrow(X)) {
      stop(
        sprintf(
          "For matrix t, nrow(t) must equal nrow(data): %d versus %d.",
          nrow(t),
          nrow(X)
        )
      )
    }

    mu_mat <- matrix(
      mu_c,
      nrow = length(mu_c),
      ncol = ncol(t)
    )

    t_safe <- pmax(t, .Machine$double.xmin)

    pr <- 1 - pnorm(
      (log(t_safe) - mu_mat) / sigma_c
    )

    pr[t <= 0] <- 1

  } else {

    if (length(t) == 1) {
      t <- rep(t, length(mu_c))
    }

    if (length(t) != length(mu_c)) {
      stop(
        sprintf(
          "Length of t must equal nrow(data): %d versus %d.",
          length(t),
          length(mu_c)
        )
      )
    }

    t_safe <- pmax(t, .Machine$double.xmin)

    pr <- 1 - pnorm(
      (log(t_safe) - mu_c) / sigma_c
    )

    pr[t <= 0] <- 1
  }

  pr[!is.finite(pr)] <- 0

  pmin(
    pmax(pr, 0),
    1
  )
}


step_survival <- function(time_grid, survival_values, t) {
  t <- as.numeric(t)

  if (length(time_grid) == 0L) {
    return(rep(1, length(t)))
  }

  idx <- findInterval(t, time_grid)
  out <- rep(1, length(t))
  use <- idx > 0L

  out[use] <- survival_values[
    pmin(idx[use], length(survival_values))
  ]
  out[t <= 0] <- 1

  pmin(pmax(out, 0), 1)
}


rowwise_step_survival <- function(time_grid, survival_matrix, t) {
  survival_matrix <- as.matrix(survival_matrix)
  n <- nrow(survival_matrix)

  if (is.matrix(t)) {
    if (nrow(t) != n) {
      stop("For matrix t, nrow(t) must equal nrow(survival_matrix).")
    }

    k <- ncol(t)
    flat_t <- as.vector(t)
    idx <- findInterval(flat_t, time_grid)
    row_idx <- rep(seq_len(n), times = k)
    out <- rep(1, length(flat_t))
    use <- idx > 0L

    out[use] <- survival_matrix[
      cbind(
        row_idx[use],
        pmin(idx[use], ncol(survival_matrix))
      )
    ]
    out[flat_t <= 0] <- 1

    return(
      matrix(
        pmin(pmax(out, 0), 1),
        nrow = n,
        ncol = k
      )
    )
  }

  if (length(t) == 1L) {
    t <- rep(t, n)
  }
  if (length(t) != n) {
    stop("Length of t must equal nrow(survival_matrix).")
  }

  idx <- findInterval(t, time_grid)
  out <- rep(1, n)
  use <- idx > 0L
  out[use] <- survival_matrix[
    cbind(
      which(use),
      pmin(idx[use], ncol(survival_matrix))
    )
  ]
  out[t <= 0] <- 1

  pmin(pmax(out, 0), 1)
}


fit_sc_model <- function(method,
                         data_fit,
                         xnames,
                         setting = NULL,
                         ntree = 1000,
                         gamma = 1.0) {
  method <- match.arg(
    method,
    c("oracle", "power_oracle", "rsf", "aft_lognormal", "km_x1", "km")
  )

  dat <- as.data.frame(data_fit)
  dat$cens_event <- 1L - as.integer(dat$event)

  if (method == "oracle") {
    if (is.null(setting)) {
      stop("setting is required for method = 'oracle'.")
    }
    return(make_oracle_sc_model(setting))
  }

  if (method == "power_oracle") {
    if (is.null(setting)) {
      stop("setting is required for method = 'power_oracle'.")
    }
    return(make_power_oracle_sc_model(setting, gamma = gamma))
  }

  if (method == "aft_lognormal") {
    fmla <- as.formula(
      paste(
        "Surv(censored_T, cens_event) ~",
        paste(xnames, collapse = " + ")
      )
    )
    fit <- survival::survreg(fmla, data = dat, dist = "lognormal")

    return(
      structure(
        list(fit = fit, xnames = xnames),
        class = "sc_aft_lognormal"
      )
    )
  }

  if (method == "km") {
    fit <- survival::survfit(
      survival::Surv(censored_T, cens_event) ~ 1,
      data = dat
    )
    return(structure(list(fit = fit), class = "sc_km"))
  }

  if (method == "km_x1") {
    if (!("X1" %in% colnames(dat))) {
      stop("km_x1 requires X1.")
    }

    fits <- setNames(vector("list", 2), c("0", "1"))
    for (g in c(0, 1)) {
      dat_g <- dat[dat$X1 == g, , drop = FALSE]
      fits[[as.character(g)]] <- if (nrow(dat_g) == 0L) {
        NULL
      } else {
        survival::survfit(
          survival::Surv(censored_T, cens_event) ~ 1,
          data = dat_g
        )
      }
    }

    return(structure(list(fits = fits), class = "sc_km_x1"))
  }

  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    stop("Install randomForestSRC before using method = 'rsf'.")
  }

  fmla <- as.formula(
    paste(
      "Surv(censored_T, cens_event) ~",
      paste(xnames, collapse = " + ")
    )
  )
  fit_dat <- dat[
    ,
    c("censored_T", "cens_event", xnames),
    drop = FALSE
  ]
  fit <- randomForestSRC::rfsrc(
    formula = fmla,
    data = fit_dat,
    ntree = ntree,
    nodesize = 15,
    nsplit = 10,
    importance = "none",
    block.size = 10
  )

  structure(
    list(fit = fit, xnames = xnames),
    class = "sc_rsf"
  )
}


sc_prob_aft_lognormal <- function(mdl0, data, t) {
  X <- as.data.frame(data)
  newdata <- X[, mdl0$xnames, drop = FALSE]
  mu <- as.numeric(predict(mdl0$fit, newdata = newdata, type = "lp"))
  sigma <- as.numeric(mdl0$fit$scale)
  n <- nrow(newdata)

  if (is.matrix(t)) {
    if (nrow(t) != n) {
      stop("For matrix t, nrow(t) must equal nrow(data).")
    }
    mu_mat <- matrix(mu, nrow = n, ncol = ncol(t))
    t_safe <- pmax(t, .Machine$double.xmin)
    out <- 1 - pnorm((log(t_safe) - mu_mat) / sigma)
    out[t <= 0] <- 1
  } else {
    if (length(t) == 1L) {
      t <- rep(t, n)
    }
    if (length(t) != n) {
      stop("Length of t must equal nrow(data).")
    }
    t_safe <- pmax(t, .Machine$double.xmin)
    out <- 1 - pnorm((log(t_safe) - mu) / sigma)
    out[t <= 0] <- 1
  }

  out[!is.finite(out)] <- 0
  pmin(pmax(out, 0), 1)
}


sc_prob_km <- function(mdl0, data, t) {
  fit <- mdl0$fit

  if (is.matrix(t)) {
    out <- step_survival(fit$time, fit$surv, as.vector(t))
    return(matrix(out, nrow = nrow(t), ncol = ncol(t)))
  }
  if (length(t) == 1L) {
    t <- rep(t, nrow(data))
  }
  step_survival(fit$time, fit$surv, t)
}


sc_prob_km_x1 <- function(mdl0, data, t) {
  X <- as.data.frame(data)
  if (!("X1" %in% colnames(X))) {
    stop("sc_km_x1 prediction requires X1.")
  }

  n <- nrow(X)
  if (is.matrix(t)) {
    if (nrow(t) != n) {
      stop("For matrix t, nrow(t) must equal nrow(data).")
    }
    out <- matrix(NA_real_, nrow = n, ncol = ncol(t))
    for (g in c(0, 1)) {
      rows <- which(X$X1 == g)
      if (length(rows) == 0L) next
      fit <- mdl0$fits[[as.character(g)]]
      out[rows, ] <- if (is.null(fit)) {
        1
      } else {
        matrix(
          step_survival(
            fit$time,
            fit$surv,
            as.vector(t[rows, , drop = FALSE])
          ),
          nrow = length(rows),
          ncol = ncol(t)
        )
      }
    }
    return(out)
  }

  if (length(t) == 1L) {
    t <- rep(t, n)
  }
  if (length(t) != n) {
    stop("Length of t must equal nrow(data).")
  }

  out <- rep(NA_real_, n)
  for (g in c(0, 1)) {
    rows <- which(X$X1 == g)
    if (length(rows) == 0L) next
    fit <- mdl0$fits[[as.character(g)]]
    out[rows] <- if (is.null(fit)) {
      1
    } else {
      step_survival(fit$time, fit$surv, t[rows])
    }
  }
  out
}


sc_prob_rsf <- function(mdl0, data, t) {
  X <- as.data.frame(data)
  newdata <- X[, mdl0$xnames, drop = FALSE]
  pred <- predict(mdl0$fit, newdata = newdata)

  rowwise_step_survival(
    time_grid = pred$time.interest,
    survival_matrix = as.matrix(pred$survival),
    t = t
  )
}


is_sc_model <- function(mdl0) {
  any(vapply(
    c("oracle_sc", "power_oracle_sc", "sc_aft_lognormal", "sc_km", "sc_km_x1", "sc_rsf"),
    function(class_name) inherits(mdl0, class_name),
    logical(1)
  ))
}


sc_prob <- function(mdl0, data, xnames, t) {
  if (inherits(mdl0, "oracle_sc")) {
    return(oracle_sc_prob(mdl0, data, t))
  }

  if (inherits(mdl0, "power_oracle_sc")) {

    oracle_mdl <- make_oracle_sc_model(
      setting = mdl0$setting
    )

    true_G <- oracle_sc_prob(
      oracle_mdl = oracle_mdl,
      data = data,
      t = t
    )

    wrong_G <- true_G ^ mdl0$gamma

    return(
      pmin(
        pmax(wrong_G, 0),
        1
      )
    )
  }

  if (inherits(mdl0, "sc_aft_lognormal")) {
    return(sc_prob_aft_lognormal(mdl0, data, t))
  }
  if (inherits(mdl0, "sc_km")) {
    return(sc_prob_km(mdl0, data, t))
  }
  if (inherits(mdl0, "sc_km_x1")) {
    return(sc_prob_km_x1(mdl0, data, t))
  }
  if (inherits(mdl0, "sc_rsf")) {
    return(sc_prob_rsf(mdl0, data, t))
  }

  Xmat <- as.matrix(data[, xnames, drop = FALSE])

  gpr_mean <- mdl0$predict(Xmat)
  gpr_sd <- mdl0$predict(Xmat, se.fit = TRUE)$se

  if (is.matrix(t)) {
    mean_mat <- matrix(gpr_mean, nrow = length(gpr_mean), ncol = ncol(t))
    sd_mat <- matrix(gpr_sd, nrow = length(gpr_sd), ncol = ncol(t))
    pr <- pnorm((-t - mean_mat) / sd_mat)
  } else {
    if (length(t) == 1) {
      t <- rep(t, length(gpr_mean))
    }

    pr <- pnorm((-t - gpr_mean) / gpr_sd)
  }

  pr
}
