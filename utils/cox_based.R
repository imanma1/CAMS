############################################
## Lower prediction bound based on Cox model
## with adaptive cutoffs
############################################

cox_based <- function(x, p, len_x, xnames,
                      data_fit, data_calib,
                      mdl0, alpha,
                      use_oracle_sc = FALSE) {

  n <- nrow(data_calib)

  ########################################
  ## Fit the survival model
  ########################################

  newdata <- data.frame(x)
  colnames(newdata) <- xnames

  fmla <- as.formula(
    paste(
      "Surv(censored_T, event) ~",
      paste(xnames, collapse = "+")
    )
  )

  fit_nonconverged <- FALSE
  fit_warning <- NULL
  fit_error <- NULL

  mdl <- tryCatch(
    withCallingHandlers(
      survreg(
        fmla,
        data = data_fit,
        dist = "weibull",
        control = survreg.control(maxiter = 1000)
      ),
      warning = function(w) {

        msg <- conditionMessage(w)

        if (grepl(
          "did not converge|Ran out of iterations",
          msg,
          ignore.case = TRUE
        )) {

          fit_nonconverged <<- TRUE
          fit_warning <<- msg

          # Suppress the original warning because we issue
          # a clearer method-level warning below.
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) {

      fit_error <<- conditionMessage(e)

      NULL
    }
  )

  ########################################
  ## Check whether the fit is usable
  ########################################

  fit_failed <- (
    is.null(mdl) ||
    fit_nonconverged ||
    !is.null(fit_error)
  )

  if (!fit_failed) {

    fit_failed <- (
      length(coef(mdl)) == 0L ||
      any(!is.finite(coef(mdl))) ||
      length(mdl$scale) != 1L ||
      !is.finite(mdl$scale) ||
      mdl$scale <= 0
    )
  }

  ########################################
  ## Optional debugging output
  ########################################

  if (isTRUE(getOption("debug_dft_fixed", FALSE))) {

    cat("\n========== DFT-ADAPTIVE SURVREG DEBUG ==========\n")
    cat("fit n:", nrow(data_fit), "\n")

    mm <- tryCatch(
      model.matrix(fmla, data = data_fit),
      error = function(e) NULL
    )

    if (!is.null(mm)) {

      cat(
        "number of model-matrix columns:",
        ncol(mm),
        "\n"
      )

      cat(
        "model-matrix rank:",
        qr(mm)$rank,
        "\n"
      )

      cat(
        "condition number:",
        kappa(mm),
        "\n"
      )
    }

    if (!is.null(mdl)) {

      cat(
        "survreg iterations:",
        mdl$iter,
        "\n"
      )

      cat(
        "non-finite coefficients:",
        sum(!is.finite(coef(mdl))),
        "/",
        length(coef(mdl)),
        "\n"
      )

      cat(
        "scale:",
        mdl$scale,
        "\n"
      )
    } else {

      cat("survreg model: NULL\n")
    }

    cat(
      "nonconvergence warning:",
      fit_nonconverged,
      "\n"
    )

    if (!is.null(fit_warning)) {
      cat(
        "warning text:",
        fit_warning,
        "\n"
      )
    }

    if (!is.null(fit_error)) {
      cat(
        "error text:",
        fit_error,
        "\n"
      )
    }

    cat(
      "fit failed:",
      fit_failed,
      "\n"
    )

    cat("=================================================\n\n")
  }

  ########################################
  ## Stop DFT-adaptive if survival fit failed
  ########################################

  if (fit_failed) {

    failure_message <- if (!is.null(fit_warning)) {

      fit_warning

    } else if (!is.null(fit_error)) {

      fit_error

    } else {

      "Non-finite or otherwise invalid survreg fit."
    }

    warning(
      paste0(
        "cox_based(): Weibull survival model failed. ",
        "Returning NA predictions for DFT-adaptive-T ",
        "and DFT-adaptive-CT. Reason: ",
        failure_message
      )
    )

    return(
      list(
        lower_bnd_qtg = rep(NA_real_, len_x),
        lower_bnd_qtl = rep(NA_real_, len_x),
        lower_bnd_qctg = rep(NA_real_, len_x),
        lower_bnd_qctl = rep(NA_real_, len_x),
        time_qt = 0,
        time_qct = 0,
        fit_failed = TRUE,
        fit_warning = failure_message
      )
    )
  }

  ########################################
  ## Cutoff for quantile of C
  ## to bound weights
  ########################################

  cens_rt <- 1 - 1 / log(n)

  ########################################
  ## DFT-adaptive-T
  ########################################

  start_time <- proc.time()[3]

  qt_res <- alpha_qt(
    mdl,
    newdata,
    data_fit,
    data_calib,
    xnames,
    alpha,
    len_x,
    cens_rt = cens_rt,
    mdl0 = mdl0,
    use_oracle_sc = use_oracle_sc
  )

  lower_bnd_qtg <- qt_res$lower_bnd_g
  lower_bnd_qtl <- qt_res$lower_bnd_l

  end_time <- proc.time()[3]
  time_qt <- end_time - start_time

  ########################################
  ## DFT-adaptive-CT
  ########################################

  start_time <- proc.time()[3]

  fit_X <- data_fit[
    ,
    xnames,
    drop = FALSE
  ]

  qc_mdl <- quantile_forest(
    fit_X,
    as.vector(data_fit$C)
  )

  qct_res <- alpha_qct(
    mdl,
    qc_mdl,
    newdata,
    data_fit,
    data_calib,
    xnames,
    alpha,
    len_x,
    cens_rt = cens_rt,
    mdl0 = mdl0,
    use_oracle_sc = use_oracle_sc
  )

  lower_bnd_qctg <- qct_res$lower_bnd_g
  lower_bnd_qctl <- qct_res$lower_bnd_l

  end_time <- proc.time()[3]
  time_qct <- end_time - start_time

  ########################################
  ## Return results
  ########################################

  return(
    list(
      lower_bnd_qtg = lower_bnd_qtg,
      lower_bnd_qtl = lower_bnd_qtl,
      lower_bnd_qctg = lower_bnd_qctg,
      lower_bnd_qctl = lower_bnd_qctl,
      time_qt = time_qt,
      time_qct = time_qct,
      fit_failed = FALSE,
      fit_warning = NULL
    )
  )
}