############################################
## Lower prediction bound based on Cox model
## with adaptive cutoffs
############################################

cox_based <- function(x, p, len_x, xnames,
                      data_fit, data_calib,
                      mdl0, alpha,
                      use_oracle_sc = FALSE) {
  n <- nrow(data_calib)
  ## Fit the survival model
  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  mdl <- survreg(fmla, data = data_fit, dist = "weibull")

  # cutoff for quantile of C
  # to bound weights
  cens_rt <- 1 - 1 / log(n)

  start_time = proc.time()[3]
  qt_res <- alpha_qt(mdl, newdata,
                     data_fit, data_calib,
                     xnames, alpha, len_x,
                     cens_rt = cens_rt,
                     mdl0 = mdl0,
                     use_oracle_sc = use_oracle_sc)
  lower_bnd_qtg <- qt_res$lower_bnd_g
  lower_bnd_qtl <- qt_res$lower_bnd_l
  end_time <- proc.time()[3]
  time_qt <- end_time - start_time

  start_time = proc.time()[3]
  ## Fit the model for C with quantile_forest (now only supports 1d)
  fit_X <- data_fit[, xnames, drop = FALSE]
  qc_mdl <- quantile_forest(fit_X, as.vector(data_fit$C))

  qct_res <- alpha_qct(mdl, qc_mdl, newdata,
                       data_fit, data_calib,
                       xnames, alpha, len_x,
                       cens_rt = cens_rt,
                       mdl0 = mdl0,
                       use_oracle_sc = use_oracle_sc)
  lower_bnd_qctg <- qct_res$lower_bnd_g
  lower_bnd_qctl <- qct_res$lower_bnd_l
  end_time <- proc.time()[3]
  time_qct <- end_time - start_time

  return(list(lower_bnd_qtg =  lower_bnd_qtg,
              lower_bnd_qtl = lower_bnd_qtl,
              lower_bnd_qctg = lower_bnd_qctg,
              lower_bnd_qctl = lower_bnd_qctl,
              time_qt = time_qt,
              time_qct = time_qct))
}