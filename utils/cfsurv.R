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
                   mdl0,
                   use_oracle_sc = FALSE) {
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

    res <- cox_censoring_prob(
      mdl0 = mdl0,
      calib = data_calib,
      test = newdata,
      xnames = xnames,
      c = c,
      ftol = ftol,
      tol = tol,
      use_oracle_sc = use_oracle_sc
    )

    pr_calib <- res$pr_calib
    pr_new <- res$pr_new

  } else {
    res <- selection_c(
      data = data_fit,
      p = p,
      n = nrow(data_fit),
      xnames = xnames,
      c_ref = c_list,
      weight_ref = NULL,
      model = model,
      type = type,
      dist = dist,
      mdl0 = mdl0,
      alpha = alpha,
      use_oracle_sc = use_oracle_sc
    )

    c <- res$c_opt

    res <- cox_censoring_prob(
      mdl0 = mdl0,
      calib = data_calib,
      test = newdata,
      xnames = xnames,
      c = c,
      ftol = ftol,
      tol = tol,
      use_oracle_sc = use_oracle_sc
    )

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

  return(list(res = res, c = c))
}


cox_censoring_prob <- function(mdl0, calib, test = NULL,
                               xnames, c,
                               ftol = .1, tol = .1,
                               use_oracle_sc = FALSE) {
  p <- length(xnames)

  if (use_oracle_sc || inherits(mdl0, "oracle_sc")) {

    pr_calib <- sc_prob(
      mdl0 = mdl0,
      data = calib,
      xnames = xnames,
      t = c
    )

    if (!is.null(test)) {
      newdata <- data.frame(test)

      if (ncol(newdata) == length(xnames)) {
        colnames(newdata) <- xnames
      }

      # In subgroup mode, xnames may exclude X1, but oracle S_C may need X1.
      # Since calib is subgroup-specific, X1 is constant and can be recovered.
      if (!("X1" %in% colnames(newdata)) && ("X1" %in% colnames(calib))) {
        x1_vals <- unique(calib$X1)

        if (length(x1_vals) == 1) {
          newdata$X1 <- x1_vals
        } else {
          stop("Oracle S_C needs X1 for test data, but X1 cannot be inferred from calibration data.")
        }
      }

      pr_new <- sc_prob(
        mdl0 = mdl0,
        data = newdata,
        xnames = xnames,
        t = c
      )

    } else {
      pr_new <- NULL
    }

    return(list(pr_calib = pr_calib, pr_new = pr_new))
  }

  ## Computing the censoring scores for the calibration data
  mean_calib <- mdl0$predict(as.matrix(calib[, xnames, drop = FALSE]))
  sd_calib <- mdl0$predict(as.matrix(calib[, xnames, drop = FALSE]), se.fit = TRUE)$se

  pr_calib <- pnorm((-c - mean_calib) / sd_calib)

  ## Computing the censoring scores for the test data
  if (!is.null(test)) {
    newdata <- data.frame(test)
    colnames(newdata) <- xnames

    mean_new <- mdl0$predict(as.matrix(newdata[, xnames, drop = FALSE]))
    sd_new <- mdl0$predict(as.matrix(newdata[, xnames, drop = FALSE]), se.fit = TRUE)$se

    pr_new <- pnorm((-c - mean_new) / sd_new)
  } else {
    pr_new <- NULL
  }

  return(list(pr_calib = pr_calib, pr_new = pr_new))
}