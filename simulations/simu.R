simu <- function(seed, setting, only_cams = FALSE,
                 n_train, n_calib, n_test,
                 xmin, xmax, alpha,
                 bernoulli_prob = 0.1,
                 use_oracle_sc = FALSE,
                 sc_method = NULL,
                 augmentation_method = c("same", "correct", "wrong", "oracle_event"),
                 homoscedastic_event = FALSE,
                 sc_ntree = 1000,
                 gamma = 1.0,
                 use_intersectional_R = FALSE) {
  set.seed(seed)
  mod <- "cox"

  if (is.null(sc_method)) {
    sc_method <- if (isTRUE(use_oracle_sc)) "oracle" else "aft_lognormal"
  }
  sc_method <- match.arg(
    sc_method,
    c("oracle", "rsf", "aft_lognormal", "km_x1", "km", "power_oracle")
  )
  augmentation_method <- match.arg(augmentation_method)
  use_oracle_sc <- identical(sc_method, "oracle")

  n_threads <- 1
  if (.Platform$OS.type != "windows") {
    library(parallel)
    library(RhpcBLASctl)
    slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    n_threads <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
    # Throttle down to 1 thread to protect mclapply from deadlocking
    blas_set_num_threads(n_threads)
    omp_set_num_threads(n_threads)
  }

  ## Initialization

  ## Generate data according to the setting
  data_obj <- model_generating_fun(n_train, n_calib, n_test,
                                   setting, xmin, xmax, bernoulli_prob,
                                   homoscedastic_event = homoscedastic_event)

  p <- data_obj$p
  xnames <- paste0("X", 1:p)
  xnames_sub <- paste0("X", 2:p) # Feature names without X1 for subgroup training

  data_fit <- data_obj$data_fit
  data_calib <- data_obj$data_calib
  data_test <- data_obj$data_test
  T_test <- data_obj$T_test

  cat(sprintf("Fitting censoring model: %s\n", sc_method))
  start_time_sc <- proc.time()[3]
  mdl0 <- fit_sc_model(
    method = sc_method,
    data_fit = data_fit,
    xnames = xnames,
    setting = setting,
    ntree = sc_ntree,
    gamma = gamma
  )
  time_mdl0 <- proc.time()[3] - start_time_sc
  cat(sprintf("Censoring model fitted in %.2f seconds.\n", time_mdl0))

  make_R_label <- function(data, use_intersectional_R = FALSE) {

  if (!all(data$X1 %in% c(0, 1))) {
    stop("X1 must contain only 0 and 1.")
  }

  if (!use_intersectional_R) {
    return(
      ifelse(
        data$X1 == 0,
        "X1_0",
        "X1_1"
      )
    )
  }

  if (!"X2" %in% colnames(data)) {
    stop(
      "X2 is required when use_intersectional_R = TRUE."
    )
  }

  paste0(
    "X1_", data$X1,
    "__X2_",
    ifelse(
      data$X2 > 0,
      "gt0",
      "le0"
    )
  )
}


if (use_intersectional_R) {

  R_levels <- c(
    "X1_0__X2_le0",
    "X1_0__X2_gt0",
    "X1_1__X2_le0",
    "X1_1__X2_gt0"
  )

} else {

  R_levels <- c(
    "X1_0",
    "X1_1"
  )
}


R_fit <- make_R_label(
  data_fit,
  use_intersectional_R
)

R_calib <- make_R_label(
  data_calib,
  use_intersectional_R
)

R_test <- make_R_label(
  data_test,
  use_intersectional_R
)

idx_test_0 <- data_test$X1 == 0
idx_test_1 <- data_test$X1 == 1

  ########################################
  ## Core Pipeline Helper Function
  ########################################
  # This function trains all 6 methods and returns their lower bounds and times
  run_pipeline <- function(sub_fit, sub_calib, sub_test,
                           xnames_to_use, alpha, seed,
                           mod, non_cams_mode = "joint") {
    # Utility function to safely extract quantiles from a coxph object
    p_sub <- length(xnames_to_use)
    len_test <- nrow(sub_test)

    x_test <- sub_test[, xnames_to_use, drop = FALSE]
    sub_data <- rbind(sub_fit, sub_calib)
    Xtrain <- sub_data[, xnames_to_use, drop = FALSE]

    if (!only_cams) {
      # 2. cfsurv_c (qt and qct)
      cat("Training cfsurv_c...\n")
      start_time <- proc.time()[3]
      lb_res <- cox_based(
        x = x_test,
        p = p_sub,
        len_x = len_test,
        xnames = xnames_to_use,
        data_fit = sub_fit,
        data_calib = sub_calib,
        mdl0 = mdl0,
        alpha = alpha,
        use_oracle_sc = use_oracle_sc
      )
      time_q <- proc.time()[3] - start_time
      cat(sprintf("cfsurv_c trained in %.2f seconds.\n", time_q))

      res1 <- lb_res$lower_bnd_qtl
      res2 <- lb_res$lower_bnd_qctl
      time_qt <- lb_res$time_qt
      time_qct <- lb_res$time_qct
      time_q_base <- time_q - (time_qt + time_qct)
      time_qt <- time_qt + time_q_base
      time_qct <- time_qct + time_q_base

      times <- c(time_qt, time_qct)
      output <- data.frame(
        "DFT-adaptive-T" = res1,
        "DFT-adaptive-CT" = res2,
        check.names = FALSE
      )

      # 3. cfsurv (qc0)
      cat("Training cfsurv (qc0)...\n")
      start_time <- proc.time()[3]
      res0 <- cfsurv(
        x = x_test,
        p = p_sub,
        len_x = len_test,
        xnames = xnames_to_use,
        data_fit = sub_fit,
        data_calib = sub_calib,
        n = nrow(sub_data),
        alpha = alpha,
        type = "quantile",
        model = mod,
        dist = "weibull",
        c_list = NULL,
        pr_list = NULL,
        pr_new_list = NULL,
        ftol = 0.1,
        tol = 0.1,
        n_tree = 100,
        mdl0 = mdl0,
        use_oracle_sc = use_oracle_sc
      )
      raw_time_qc0 <- proc.time()[3] - start_time
      cat(sprintf("cfsurv (qc0) trained in %.2f seconds.\n", raw_time_qc0))
      time_qc0 <- raw_time_qc0
      times <- c(times, time_qc0)
      if (length(res0$res) == 0) {
        cat("  -> [WARNING] cfsurv returned empty predictions. Filling with NAs.\n")
        output[["DFT-fixed"]] <- rep(NA, nrow(output))
      } else {
        output[["DFT-fixed"]] <- res0$res
      }

      # ==========================================================
      # 4. Vanilla CQR on the original time scale
      # ==========================================================
      cat("Training vanilla CQR on original scale...\n")

      start_time <- proc.time()[3]

      res_raw <- lapply(
        alpha,
        cqr,
        x = x_test,
        Xtrain = Xtrain,
        Ytrain = sub_data$censored_T,
        I_fit = seq_len(nrow(sub_fit)),
        seed = seed + 7
      )

      res_raw <- do.call(
        rbind,
        lapply(res_raw, as.data.frame)
      )

      output[["Vanilla CQR"]] <- as.numeric(res_raw[, 1])

      vanilla_cqr_time <- proc.time()[3] - start_time

      times <- c(
        times,
        vanilla_cqr_time
      )

      cat(
        sprintf(
          "Vanilla CQR on original scale trained in %.2f seconds.\n",
          vanilla_cqr_time
        )
      )


      # ==========================================================
      # 4b. Vanilla CQR on the log-time scale
      # ==========================================================
      cat("Training vanilla CQR on log scale...\n")

      start_time <- proc.time()[3]

      log_censored_time <- log(
        pmax(
          sub_data$censored_T,
          .Machine$double.xmin
        )
      )

      res_log <- lapply(
        alpha,
        cqr,
        x = x_test,
        Xtrain = Xtrain,
        Ytrain = log_censored_time,
        I_fit = seq_len(nrow(sub_fit)),
        seed = seed + 17
      )

      res_log <- do.call(
        rbind,
        lapply(res_log, as.data.frame)
      )

      log_lower_bound <- as.numeric(res_log[, 1])

      # Avoid numerical overflow before exponentiating
      log_lower_bound <- pmin(
        log_lower_bound,
        log(.Machine$double.xmax)
      )

      log_lower_bound <- pmax(
        log_lower_bound,
        log(.Machine$double.xmin)
      )

      output[["Vanilla CQR-log"]] <- exp(log_lower_bound)

      vanilla_cqr_log_time <- proc.time()[3] - start_time

      times <- c(
        times,
        vanilla_cqr_log_time
      )

      cat(
        sprintf(
          "Vanilla CQR on log scale trained in %.2f seconds.\n",
          vanilla_cqr_log_time
        )
      )

      # 5. Cox Model
      # cat("Training Cox Model...\n")
      # start_time <- proc.time()[3]
      # fmla <- as.formula(paste("Surv(censored_T, event) ~", paste(xnames_to_use, collapse="+")))
      # mdl <- coxph(fmla, data = sub_data)
      # # mclapply distributes the rows across all available cores
      # if (.Platform$OS.type != "windows") {
      #   original_threads <- blas_get_num_procs()
      #   # Throttle down to 1 thread to protect mclapply from deadlocking
      #   blas_set_num_threads(1)
      #   omp_set_num_threads(1)
      # }
      # cox_res_list <- mclapply(1:nrow(sub_test), function(i) {
      #   extract_quant(mdl, sub_test[i, , drop=FALSE], alpha)
      # }, mc.cores = n_threads)

      # if (.Platform$OS.type != "windows") {
      #   # Restore the original number of threads after parallel processing
      #   blas_set_num_threads(original_threads)
      #   omp_set_num_threads(original_threads)
      # }
      # cox_res <- unlist(cox_res_list)
      # output$cox.bnd <- cox_res
      # cox_time <- proc.time()[3] - start_time
      # times <- c(times, cox_time)
      # cat(sprintf("Cox Model trained in %.2f seconds.\n", cox_time))
      
      # 6. Random Forest

      # extract_quant <- function(mdl, x, alpha) {
      # res <- summary(survfit(mdl, newdata = x))
      # time_point <- res$time
      # survcdf <- 1 - res$surv
      # if(sum(survcdf >= alpha)==0){
      #   quant = max(time_point)
      # }else{
      #   quant <- time_point[min(which(survcdf >= alpha))]
      # }
      # return(quant)
      # }

      # cat("Training Random Forest...\n")
      # start_time <- proc.time()[3]
      # ntree <- 1000
      # nodesize <- 80
      # fmla_rf <- as.formula(paste("censored_T ~", paste(xnames_to_use, collapse="+")))
      # mdl <- crf.km(fmla_rf, ntree = ntree, nodesize = nodesize,
      #               data_train = sub_data[, c(xnames_to_use, "censored_T", "event"), drop=FALSE], 
      #               data_test = sub_test[, xnames_to_use, drop=FALSE], 
      #               yname = 'censored_T', iname = 'event', tau = alpha, method = "grf")
      # output$rf.bnd <- mdl$predicted
      # rf_time <- proc.time()[3] - start_time
      # times <- c(times, rf_time)
      # cat(sprintf("Random Forest trained in %.2f seconds.\n", rf_time))
    }
    else {
      times <- NA
      output <- NA
    }

    return(list(output = output, times = times, mdl0 = mdl0))
  }

  compute_metrics <- function(
      output_df,
      times_vec = NULL,
      suffix_label = NULL,
      calibration_diagnostics = NULL,
      use_intersectional_R = FALSE
  ) {

    output_df[] <- lapply(
      output_df,
      function(z) {
        z[is.finite(z)] <- pmax(
          z[is.finite(z)],
          0
        )
        z
      }
    )

    method_names_original <- colnames(output_df)
    method_names <- method_names_original

    if (!is.null(suffix_label)) {
      method_names <- paste(
        method_names,
        suffix_label
      )
    }

    if (is.null(times_vec)) {
      times_vec <- rep(
        NA_real_,
        ncol(output_df)
      )
    }

    if (length(times_vec) == 1L) {
      times_vec <- rep(
        times_vec,
        ncol(output_df)
      )
    }

    if (length(times_vec) != ncol(output_df)) {
      stop(
        "times_vec must have length 1 or ncol(output_df)."
      )
    }

    safe_group_mean <- function(x, index) {

      if (!any(index)) {
        return(NA_real_)
      }

      mean(
        x[index],
        na.rm = TRUE
      )
    }

    weighted_finite_mean <- function(
        values,
        weights
    ) {

      values <- suppressWarnings(
        as.numeric(values)
      )

      weights <- suppressWarnings(
        as.numeric(weights)
      )

      valid <- is.finite(values) &
        is.finite(weights) &
        weights > 0

      if (!any(valid)) {
        return(NA_real_)
      }

      weighted.mean(
        values[valid],
        weights[valid]
      )
    }

    # ----------------------------------------------------------
    # Define the evaluation groups
    # ----------------------------------------------------------

    if (use_intersectional_R) {

      group_specs <- list(
        list(
          key = "X1_0__X2_le0",
          label = "x_1 = 0 and x_2 <= 0",
          index = data_test$X1 == 0 &
            data_test$X2 <= 0
        ),
        list(
          key = "X1_0__X2_gt0",
          label = "x_1 = 0 and x_2 > 0",
          index = data_test$X1 == 0 &
            data_test$X2 > 0
        ),
        list(
          key = "X1_1__X2_le0",
          label = "x_1 = 1 and x_2 <= 0",
          index = data_test$X1 == 1 &
            data_test$X2 <= 0
        ),
        list(
          key = "X1_1__X2_gt0",
          label = "x_1 = 1 and x_2 > 0",
          index = data_test$X1 == 1 &
            data_test$X2 > 0
        )
      )

      R_scheme <- "X1-by-X2-sign"

    } else {

      group_specs <- list(
        list(
          key = "X1_0",
          label = "x_1 = 0",
          index = data_test$X1 == 0
        ),
        list(
          key = "X1_1",
          label = "x_1 = 1",
          index = data_test$X1 == 1
        )
      )

      R_scheme <- "X1"
    }

    # ----------------------------------------------------------
    # Coverage and lower-bound metrics
    # ----------------------------------------------------------

    cov_marg <- apply(
      output_df,
      2,
      function(x) {
        mean(
          T_test >= x,
          na.rm = TRUE
        )
      }
    )

    simulen <- apply(
      output_df,
      2,
      mean,
      na.rm = TRUE
    )

    group_coverage <- list()
    group_lower_bound_mean <- list()

    for (group_spec in group_specs) {

      current_index <- group_spec$index
      current_key <- group_spec$key

      group_coverage[[current_key]] <- apply(
        output_df,
        2,
        function(x) {
          safe_group_mean(
            T_test >= x,
            current_index
          )
        }
      )

      group_lower_bound_mean[[current_key]] <- apply(
        output_df,
        2,
        function(x) {
          safe_group_mean(
            x,
            current_index
          )
        }
      )
    }

    metrics <- data.frame(
      "method" = method_names,
      "setting" = setting,
      "censoring model" = sc_method,
      "augmentation model" = augmentation_method,
      "event-time scale" = if (
        homoscedastic_event
      ) {
        "constant"
      } else {
        "group-varying"
      },
      "R scheme" = R_scheme,
      "Marginal coverage" = cov_marg,
      "true test-set miscoverage" = 1 - cov_marg,
      "lower bound values mean" = simulen,
      "computation time" = times_vec,
      check.names = FALSE,
      row.names = NULL
    )

    # Add the group-specific coverage columns.
    for (group_spec in group_specs) {

      metrics[[
        paste0(
          "group coverage for ",
          group_spec$label
        )
      ]] <- group_coverage[[
        group_spec$key
      ]]
    }

    # Add the group-specific lower-bound columns.
    for (group_spec in group_specs) {

      metrics[[
        paste0(
          "lower bound mean for ",
          group_spec$label
        )
      ]] <- group_lower_bound_mean[[
        group_spec$key
      ]]
    }

    # ----------------------------------------------------------
    # Calibration diagnostic columns
    # ----------------------------------------------------------

    overall_diagnostic_columns <- c(
      "selected calibration level",
      "estimated calibration risk",
      "absolute calibration-to-test risk error",
      "fraction of censoring probabilities truncated at eta"
    )

    group_diagnostic_columns <- unlist(
      lapply(
        group_specs,
        function(group_spec) {
          c(
            paste0(
              "selected calibration level for ",
              group_spec$label
            ),
            paste0(
              "estimated calibration risk for ",
              group_spec$label
            ),
            paste0(
              "absolute calibration-to-test risk error for ",
              group_spec$label
            )
          )
        }
      ),
      use.names = FALSE
    )

    metrics[
      c(
        overall_diagnostic_columns,
        group_diagnostic_columns
      )
    ] <- NA_real_

    # ----------------------------------------------------------
    # Fill the calibration diagnostics
    # ----------------------------------------------------------

    if (!is.null(calibration_diagnostics)) {

      if (
        !"method" %in%
          colnames(calibration_diagnostics)
      ) {
        stop(
          "calibration_diagnostics must contain a method column."
        )
      }

      diagnostic_group_column <- if (
        "R_label" %in%
          colnames(calibration_diagnostics)
      ) {
        "R_label"
      } else if (
        "subgroup" %in%
          colnames(calibration_diagnostics)
      ) {
        "subgroup"
      } else {
        stop(
          paste0(
            "calibration_diagnostics must contain either ",
            "R_label or subgroup."
          )
        )
      }

      diagnostic_R_labels <- as.character(
        calibration_diagnostics[[
          diagnostic_group_column
        ]]
      )

      # Backward compatibility for the original two groups.
      if (!use_intersectional_R) {

        diagnostic_R_labels[
          diagnostic_R_labels %in%
            c("0", "X1=0", "X1_0")
        ] <- "X1_0"

        diagnostic_R_labels[
          diagnostic_R_labels %in%
            c("1", "X1=1", "X1_1")
        ] <- "X1_1"
      }

      calibration_diagnostics$R_label_internal <-
        diagnostic_R_labels

      for (
        row_idx in seq_len(
          nrow(metrics)
        )
      ) {

        method_name <- method_names_original[
          row_idx
        ]

        method_diag <- calibration_diagnostics[
          calibration_diagnostics$method ==
            method_name,
          ,
          drop = FALSE
        ]

        if (nrow(method_diag) == 0L) {
          next
        }

        if (
          !"n_calib" %in%
            colnames(method_diag)
        ) {
          method_diag$n_calib <- 1
        }

        # Overall calibration summaries,
        # weighted by the calibration group sizes.
        metrics[[
          "selected calibration level"
        ]][row_idx] <- weighted_finite_mean(
          method_diag$
            selected_calibration_level,
          method_diag$n_calib
        )

        metrics[[
          "estimated calibration risk"
        ]][row_idx] <- weighted_finite_mean(
          method_diag$
            estimated_calibration_risk,
          method_diag$n_calib
        )

        metrics[[
          paste0(
            "fraction of censoring ",
            "probabilities truncated at eta"
          )
        ]][row_idx] <- weighted_finite_mean(
          method_diag$
            fraction_censoring_probabilities_truncated,
          method_diag$n_calib
        )

        metrics[[
          "absolute calibration-to-test risk error"
        ]][row_idx] <- abs(
          metrics[[
            "estimated calibration risk"
          ]][row_idx] -
            metrics[[
              "true test-set miscoverage"
            ]][row_idx]
        )

        # Group-specific diagnostics.
        for (group_spec in group_specs) {

          group_diag <- method_diag[
            method_diag$R_label_internal ==
              group_spec$key,
            ,
            drop = FALSE
          ]

          if (nrow(group_diag) == 0L) {
            next
          }

          selected_level_column <- paste0(
            "selected calibration level for ",
            group_spec$label
          )

          estimated_risk_column <- paste0(
            "estimated calibration risk for ",
            group_spec$label
          )

          risk_error_column <- paste0(
            paste0(
              "absolute calibration-to-test ",
              "risk error for "
            ),
            group_spec$label
          )

          coverage_column <- paste0(
            "group coverage for ",
            group_spec$label
          )

          metrics[[
            selected_level_column
          ]][row_idx] <- weighted_finite_mean(
            group_diag$
              selected_calibration_level,
            group_diag$n_calib
          )

          metrics[[
            estimated_risk_column
          ]][row_idx] <- weighted_finite_mean(
            group_diag$
              estimated_calibration_risk,
            group_diag$n_calib
          )

          metrics[[
            risk_error_column
          ]][row_idx] <- abs(
            metrics[[
              estimated_risk_column
            ]][row_idx] -
              (
                1 -
                  metrics[[
                    coverage_column
                  ]][row_idx]
              )
          )
        }
      }
    }

    metrics
  }

  ########################################
  ## APPROACH 1: Joint Modeling
  ########################################
  cat("========== Executing Approach 1: Joint Modeling ==========\n")
  res_joint <- run_pipeline(data_fit, data_calib, data_test,
                            xnames, alpha, seed, mod, "joint")

  ########################################
  ## APPROACH 2: Subgroup Modeling
  ########################################

  if (!only_cams) {

    subgroup_results <- vector(
      "list",
      length(R_levels)
    )

    names(subgroup_results) <- R_levels

    output_subgroup <- NULL

    for (r in R_levels) {

      idx_fit_r <- R_fit == r
      idx_calib_r <- R_calib == r
      idx_test_r <- R_test == r

      group_sizes <- c(
        fit = sum(idx_fit_r),
        calib = sum(idx_calib_r),
        test = sum(idx_test_r)
      )

      cat(
        sprintf(
          paste0(
            "\n========== Executing subgroup %s ",
            "(fit=%d, calib=%d, test=%d) ==========\n"
          ),
          r,
          group_sizes["fit"],
          group_sizes["calib"],
          group_sizes["test"]
        )
      )

      if (any(group_sizes == 0L)) {
        stop(
          sprintf(
            paste0(
              "Subgroup %s is empty in at least one ",
              "of fit/calibration/test samples."
            ),
            r
          )
        )
      }

      res_r <- run_pipeline(
        sub_fit = data_fit[
          idx_fit_r,
          ,
          drop = FALSE
        ],
        sub_calib = data_calib[
          idx_calib_r,
          ,
          drop = FALSE
        ],
        sub_test = data_test[
          idx_test_r,
          ,
          drop = FALSE
        ],
        xnames_to_use = xnames_sub,
        alpha = alpha,
        seed = seed,
        mod = mod,
        non_cams_mode = "subgroup"
      )

      subgroup_results[[r]] <- res_r

      if (is.null(output_subgroup)) {

        output_subgroup <- as.data.frame(
          matrix(
            NA_real_,
            nrow = nrow(data_test),
            ncol = ncol(res_r$output)
          ),
          check.names = FALSE
        )

        colnames(output_subgroup) <- colnames(
          res_r$output
        )

      } else if (
        !identical(
          colnames(output_subgroup),
          colnames(res_r$output)
        )
      ) {

        stop(
          sprintf(
            "Output columns are inconsistent for subgroup %s.",
            r
          )
        )
      }

      output_subgroup[
        idx_test_r,
        colnames(res_r$output)
      ] <- res_r$output
    }

    # Sum computation times across all two or four subgroup fits.
    times_subgroup <- Reduce(
      `+`,
      lapply(
        subgroup_results,
        function(result) result$times
      )
    )

    if (anyNA(output_subgroup)) {
      warning(
        "Some subgroup baseline predictions are NA."
      )
    }
  }


  ########################################
  ## APPROACH 3: CAMS
  ########################################
  cat("========== Executing Approach 3: CAMS ==========\n")
  start_time_cams <- proc.time()[3]

  cams_res <- cams(
    x = data_test[, xnames, drop = FALSE],
    p = p,
    len_x = nrow(data_test),
    xnames = xnames,
    data_fit = data_fit,
    data_calib = data_calib,
    mdl0 = res_joint$mdl0,
    alpha = alpha,
    use_oracle_sc = use_oracle_sc,
    augmentation_method = augmentation_method,
    setting = setting,
    homoscedastic_event = homoscedastic_event,
    use_intersectional_R = use_intersectional_R
  )

  # local_cams <- data.frame(check.names = FALSE)
  # for(audit_fraction in c(0.1, 0.2, 0.3, 0.4)) {
  #   calib_split <- local_stratified_calibration_split(
  #     data_calib = data_calib,
  #     audit_fraction = audit_fraction,
  #     split_seed = seed,
  #     group_name = "X1"
  #   )
  #   data_audit <- calib_split$audit
  #   data_final_calib <- calib_split$final_calibration

  #   start_time_local_cams <- proc.time()[3]

  #   local_cams_res <- cams_local_ipcw(
  #     x = data_test[, xnames, drop = FALSE],
  #     p = p,
  #     len_x = nrow(data_test),
  #     xnames = xnames,
  #     data_fit = data_fit,
  #     data_calib = data_calib,
  #     mdl0 = mdl0,
  #     alpha = alpha,
  #     data_audit = data_audit,
  #     data_local_calib = data_final_calib,

  #     calibration_split_seed = seed,

  #     mapping_log_dir = file.path(
  #       paste0("../local_mapping", audit_fraction),
  #       paste0("setting_", setting),
  #       paste0("seed_", seed)
  #     )
  #   )

  # cams_res0.2 <- cams(
  #   x = data_test[, xnames, drop = FALSE],
  #   p = p,
  #   len_x = nrow(data_test),
  #   xnames = xnames,
  #   data_fit = data_fit,
  #   data_calib = data_final_calib0.2,
  #   mdl0 = res_joint$mdl0,
  #   alpha = alpha,
  #   use_oracle_sc = use_oracle_sc,
  #   augmentation_method = augmentation_method,
  #   setting = setting,
  #   homoscedastic_event = homoscedastic_event
  # )

  #   # df_cams0.2 <- compute_metrics(
  #   #   cams_res0.2$output,
  #   #   times_vec = time_cams,
  #   #   suffix_label = "(20%)",
  #   #   calibration_diagnostics = cams_res0.2$diagnostics
  #   # )

  #   local_df_cams <- compute_metrics(
  #     local_cams_res$output,
  #     times_vec = proc.time()[3] - start_time_local_cams + time_mdl0,
  #     suffix_label = sprintf("(%d%%)", audit_fraction * 100),
  #   )

  #   local_cams <- rbind(local_cams, local_df_cams)
  # }

  time_cams <- proc.time()[3] - start_time_cams + time_mdl0
  cat(sprintf("CAMS trained in %.2f seconds.\n", time_cams))

  df_cams <- compute_metrics(
    cams_res$output,
    times_vec = time_cams,
    suffix_label = NULL,
    calibration_diagnostics = cams_res$diagnostics,
    use_intersectional_R = use_intersectional_R
  )

  ########################################
  ## Compute & Bind Final Results
  ########################################
  if (!only_cams) {
    df_joint <- compute_metrics(
      res_joint$output,
      times_vec = res_joint$times,
      suffix_label = "(Joint)",
      use_intersectional_R = use_intersectional_R
    )
    df_subgroup <- compute_metrics(
      output_subgroup,
      times_vec = times_subgroup,
      suffix_label = "(Subgroup)",
      use_intersectional_R = use_intersectional_R
    )
    # Append CAMS to the final CSV output
    simu_out <- rbind(df_cams, df_joint, df_subgroup)
  } else {
    # simu_out <- local_cams
    simu_out <- df_cams
  }

  rownames(simu_out) <- NULL
  return(simu_out)
}