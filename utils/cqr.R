#' Conformalized Quantile Regression
#'
#' @export

#' Conformalized Quantile Regression
#'
#' @export

cqr <- function(x, Xtrain, Ytrain,
                alpha = 0.05,
                I_fit = NULL,
                seed = 24601) {

  X <- Xtrain
  Y <- Ytrain

  if (is.null(dim(X)[1])) {
    n <- length(X)
    p <- 1
  } else {
    n <- dim(X)[1]
    p <- dim(X)[2]
  }

  xnames <- paste0("X", 1:p)

  data <- as.data.frame(
    cbind(Y, X)
  )

  colnames(data) <- c(
    "Y",
    xnames
  )

  ########################################
  ## Divide fitting and calibration data
  ########################################

  n_train <- n / 2
  n_calib <- n - n_train

  if (is.null(I_fit)) {
    I_fit <- sample(
      1:n,
      n_train,
      replace = FALSE
    )
  }

  data_fit <- data[
    I_fit,
    ,
    drop = FALSE
  ]

  data_calib <- data[
    -I_fit,
    ,
    drop = FALSE
  ]

  newdata <- data.frame(x)
  colnames(newdata) <- xnames

  ########################################
  ## Quantile-regression formula
  ########################################

  fmla <- as.formula(
    paste(
      "Y ~",
      paste(
        xnames,
        collapse = "+"
      )
    )
  )

  ########################################
  ## Check whether the model is estimable
  ########################################

  design_matrix <- tryCatch(
    model.matrix(
      fmla,
      data = data_fit
    ),
    error = function(e) NULL
  )

  fit_failed <- (
    is.null(design_matrix) ||
    nrow(design_matrix) < ncol(design_matrix)
  )

  if (!fit_failed) {

    fit_failed <- (
      qr(design_matrix)$rank <
        ncol(design_matrix)
    )
  }

  if (fit_failed) {

    warning(
      sprintf(
        paste0(
          "cqr(): quantile-regression design matrix ",
          "is underdetermined or rank deficient ",
          "(n_fit = %d, parameters = %d). ",
          "Returning NA predictions."
        ),
        nrow(data_fit),
        if (is.null(design_matrix)) {
          p + 1L
        } else {
          ncol(design_matrix)
        }
      )
    )

    return(
      rep(
        NA_real_,
        nrow(newdata)
      )
    )
  }

  ########################################
  ## Fit quantile regression
  ########################################

  mdl <- tryCatch(
    rq(
      fmla,
      data = data_fit,
      tau = alpha
    ),
    error = function(e) {

      warning(
        paste0(
          "cqr(): quantile-regression fit failed: ",
          conditionMessage(e),
          ". Returning NA predictions."
        )
      )

      NULL
    }
  )

  if (is.null(mdl)) {

    return(
      rep(
        NA_real_,
        nrow(newdata)
      )
    )
  }

  ########################################
  ## Calibration predictions
  ########################################

  quant <- tryCatch(
    predict(
      mdl,
      newdata = data_calib,
      type = "quantile"
    ),
    error = function(e) NULL
  )

  ########################################
  ## Test predictions
  ########################################

  new_quant <- tryCatch(
    predict(
      mdl,
      newdata = newdata,
      type = "quantile"
    ),
    error = function(e) NULL
  )

  if (
    is.null(quant) ||
    is.null(new_quant) ||
    any(!is.finite(quant)) ||
    any(!is.finite(new_quant))
  ) {

    warning(
      paste0(
        "cqr(): prediction failed or produced ",
        "non-finite values. Returning NA predictions."
      )
    )

    return(
      rep(
        NA_real_,
        nrow(newdata)
      )
    )
  }

  ########################################
  ## Conformal calibration
  ########################################

  score <- quant - data_calib$Y

  corr_term <- quantile(
    c(
      score,
      Inf
    ),
    1 - alpha
  )

  lower_bnd <- new_quant - corr_term

  return(
    as.numeric(lower_bnd)
  )
}
