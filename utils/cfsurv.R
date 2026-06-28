### added field mdl0 to the function and changed res variable to mdl0.
### changed the default value of ref_length from 100 to 30.

#' Predictive confidence interval for survival data
#'
#' The main function to generate a predictive conformal confidence interval for a unit's survival time.
#'
#' @param x a vector of the covariate for test point. 
#' @param c the censoring time for the test point.
#' @param Xtrain a n-by-p matrix of the covariate of the training data.
#' @param C a length n vector of the censoring time of the training data.
#' @param event a length n vector of indicators if the time observed is censored. TRUE corresponds to NOT censored, and FALSE censored.
#' @param time  a vevtor of length n, containing the observed survival time.
#' @param alpha a number between 0 and 1, speciifying the miscoverage rate.
#' @param seed an integer random seed (default: 24601).
#' @param model Options include "cox", "randomforest", "Powell", "Portnoy" and "PengHuang". This determines the model used to fit the condditional quantile (default: "cox").
#' @param dist either "weibull", "exponential" or "gaussian" (default: "weibull"). The distribution of T used in the cox model. 
#' @param h the bandwidth for the local confidence interval. Default is 1.
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @examples
#' # Generate data
#' n <- 500
#' X <- runif(n,0,2)
#' T <- exp(X+rnorm(n,0,1))
#' R <- rexp(n,rate = 0.01)
#' event <- T<=R
#' time <- pmin(T,R)
#' data <- data.frame(X=X,R=R,event=event,censored_T=censored_T)
#' # Prediction point
#' x <- seq(0,2,by=.4)
#' r <- 2
#' # Run cfsurv
#' res <- cfsurv(x,r,X,R,event,time,alpha=0.1,model="cox")
#'
#' @export

# function to construct conformal confidence interval
cfsurv <- function(x, p, len_x, xnames,
                   data_fit, data_calib, n,
                   alpha = 0.05,
                   type = "quantile",
                   model = "cox",
                   dist = "weibull",
                   c_list = NULL,
                   pr_list = NULL,
                   pr_new_list = NULL,
                   ftol = 0.1,
                   tol = 0.1,
                   n_tree = 100,
                   mdl0) {
  ## Check if the required packages are installed
  ## Solution found from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
  list.of.packages <- c("ggplot2",
                        "quantreg",
                        "grf",
                        "quantregForest",
                        "randomForestSRC",
                        "survival",
                        "tidyverse",
                        "fishmethods",
                        "foreach",
                        "doParallel",
                        "GauPro",
                        "gbm",
                        "np")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
  suppressPackageStartupMessages(res <- lapply(X=list.of.packages,FUN=require,character.only=TRUE))

  ## Split the data into the training set and the calibration set
  newdata <- data.frame(x)
  colnames(newdata) <- xnames

  ## If c is not specified, select c automatically
  if (!is.null(pr_list) || !is.null(pr_new_list)) {
    stop("Precomputed pr_list/pr_new_list are not supported in the refactored split-explicit version.")
  }

  if (is.null(c_list)) {
    ref_length <- 100
    c_list <- seq(min(data_fit$C), max(data_fit$C), length = ref_length)
  }

  if (length(c_list) == 1) {
    c <- c_list
    res <- cox_censoring_prob(mdl0, data_calib, newdata, xnames, c, ftol, tol)
    pr_calib <- res$pr_calib
    pr_new <- res$pr_new
  } else {
    res <- selection_c(data_fit, p, nrow(data_fit), xnames,
                       c_ref = c_list, weight_ref = NULL,
                       model = model, type = type, dist = dist,
                       mdl0 = mdl0, alpha = alpha)
    c <- res$c_opt
    res <- cox_censoring_prob(mdl0, data_calib, newdata, xnames, c, ftol, tol)
    pr_calib <- res$pr_calib
    pr_new <- res$pr_new
  }

  ## Computing the weight for the calibration data and the test data
  weight_calib <- 1 / pr_calib
  weight_new <- 1 / pr_new

  ## Run the main function and gather resutls
  res <- cox0_based(x, p, len_x, xnames,
                  c, alpha,
                  data_fit,
                  data_calib,
                  type,
                  dist,
                  weight_calib,
                  weight_new,
                  ftol,
                  tol)

  return (list(res = res, c = c))

}

cox_censoring_prob <- function(gpr_mdl, calib, test = NULL,
                               xnames, c,
                               ftol = .1, tol = .1) {
  p <- length(xnames)

  ## Computing the censoring scores for the calibration data
  mean_calib <- gpr_mdl$predict(as.matrix(calib[, xnames, drop = FALSE]))
  sd_calib <- gpr_mdl$predict(as.matrix(calib[, xnames, drop = FALSE]), se.fit = TRUE)$se

  pr_calib <- pnorm((-c - mean_calib) / sd_calib)

  ## Computing the censoring scores for the test data
  if (!is.null(test)) {
    newdata <- data.frame(test)
    colnames(newdata) <- xnames
    mean_new <- gpr_mdl$predict(as.matrix(newdata[, xnames, drop = FALSE]))
    sd_new <- gpr_mdl$predict(as.matrix(newdata[, xnames, drop = FALSE]), se.fit = TRUE)$se
    pr_new <- pnorm((-c - mean_new) / sd_new)
  } else {
    pr_new = NULL
  }
  return(list(pr_calib = pr_calib, pr_new = pr_new))
}
