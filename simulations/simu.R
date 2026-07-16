simu <- function(seed, setting, only_cams = FALSE,
                 n_train, n_calib, n_test,
                 xmin, xmax, alpha,
                 bernoulli_prob = 0.1,
                 use_oracle_sc = FALSE,
                 sc_method = NULL,
                 augmentation_method = c("same", "correct", "wrong", "oracle_event"),
                 homoscedastic_event = FALSE,
                 sc_ntree = 1000,
                 gamma = 1.0) {
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

  if (
    setting %in% c(
      "cams_vs_vanilla_lower_tail_hd_mild",
      "cams_vs_vanilla_lower_tail_hd_main",
      "cams_vs_vanilla_lower_tail_hd_strong"
    )
  ) {

    censoring_prob_at_observed_time <- sc_prob(
      mdl0 = mdl0,
      data = data_calib,
      xnames = xnames,
      t = data_calib$censored_T
    )

    cat(sprintf("\n%s censoring probability diagnostics:\n", sc_method))

    print(
      summary(
        censoring_prob_at_observed_time
      )
    )

    cat(
      sprintf(
        "Number of non-finite censoring probabilities: %d\n",
        sum(!is.finite(censoring_prob_at_observed_time))
      )
    )

    cat(
      sprintf(
        "Number outside [0,1]: %d\n",
        sum(
          censoring_prob_at_observed_time < 0 |
            censoring_prob_at_observed_time > 1,
          na.rm = TRUE
        )
      )
    )
  }

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

  idx_test_0 <- data_test$X1 == 0
  idx_test_1 <- data_test$X1 == 1

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
    cat("========== Executing Approach 2: Subgroup X1 = 0 ==========\n")
    res_0 <- run_pipeline(data_fit[data_fit$X1 == 0, , drop=FALSE],
                          data_calib[data_calib$X1 == 0, , drop=FALSE],
                          data_test[data_test$X1 == 0, , drop=FALSE],
                          xnames_sub, alpha, seed, mod, "subgroup0")
                                  
    cat("========== Executing Approach 2: Subgroup X1 = 1 ==========\n")
    res_1 <- run_pipeline(data_fit[data_fit$X1 == 1, , drop=FALSE],
                          data_calib[data_calib$X1 == 1, , drop=FALSE],
                          data_test[data_test$X1 == 1, , drop=FALSE],
                          xnames_sub, alpha, seed, mod, "subgroup1")
    
    # Merge subgroup outputs to map exactly to the data_test row order
    output_subgroup <- data.frame(matrix(ncol = ncol(res_0$output), nrow = nrow(data_test)))
    colnames(output_subgroup) <- colnames(res_0$output)
    output_subgroup[idx_test_0, ] <- res_0$output
    output_subgroup[idx_test_1, ] <- res_1$output
    times_subgroup <- res_0$times + res_1$times
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
    homoscedastic_event = homoscedastic_event
  )
  time_cams <- proc.time()[3] - start_time_cams + time_mdl0
  cat(sprintf("CAMS trained in %.2f seconds.\n", time_cams))

  # start_time_cams <- proc.time()[3]
  # local_cams_res <- new_cams(
  #   x = data_test[, xnames, drop = FALSE],
  #   p = p,
  #   len_x = nrow(data_test),
  #   xnames = xnames,
  #   data_fit = data_fit,
  #   data_calib = data_calib,
  #   mdl0 = res_joint$mdl0,
  #   alpha = alpha,
  #   use_oracle_sc = use_oracle_sc
  # )
  # local_time_cams <- proc.time()[3] - start_time_cams
  # cat(sprintf("Local CAMS trained in %.2f seconds.\n", local_time_cams))

  # save_local_cams_info_allseeds(
  #   cams_res = local_cams_res,
  #   setting = setting,
  #   seed = seed
  # )

  compute_metrics <- function(output_df, times_vec = NULL, suffix_label = NULL,
                              calibration_diagnostics = NULL) {

    output_df[] <- lapply(
      output_df,
      function(z) {
        z[is.finite(z)] <- pmax(z[is.finite(z)], 0)
        z
      }
    )

    method_names <- colnames(output_df)

    if (!is.null(suffix_label)) {
      method_names <- paste(method_names, suffix_label)
    }

    if (is.null(times_vec)) {
      times_vec <- rep(NA_real_, ncol(output_df))
    }

    if (length(times_vec) == 1) {
      times_vec <- rep(times_vec, ncol(output_df))
    }

    cov_marg <- apply(output_df, 2, function(x) {
      mean(T_test >= x, na.rm = TRUE)
    })

    cov_grp0 <- apply(output_df, 2, function(x) {
      mean(T_test[idx_test_0] >= x[idx_test_0], na.rm = TRUE)
    })

    cov_grp1 <- apply(output_df, 2, function(x) {
      mean(T_test[idx_test_1] >= x[idx_test_1], na.rm = TRUE)
    })

    simulen <- apply(output_df, 2, mean, na.rm = TRUE)

    simulen_grp0 <- apply(output_df, 2, function(x) {
      mean(x[idx_test_0], na.rm = TRUE)
    })

    simulen_grp1 <- apply(output_df, 2, function(x) {
      mean(x[idx_test_1], na.rm = TRUE)
    })

    metrics <- data.frame(
      "method"                       = method_names,
      "setting"                      = setting,
      "censoring model"              = sc_method,
      "augmentation model"           = augmentation_method,
      "event-time scale"             = if (homoscedastic_event) "constant" else "group-varying",
      "group coverage for x_1 = 0"   = cov_grp0,
      "group coverage for x_1 = 1"   = cov_grp1,
      "Marginal coverage"            = cov_marg,
      "true test-set miscoverage"    = 1 - cov_marg,
      "lower bound mean for x_1 = 0" = simulen_grp0,
      "lower bound mean for x_1 = 1" = simulen_grp1,
      "lower bound values mean"      = simulen,
      "computation time"             = times_vec,
      check.names = FALSE,
      row.names = NULL
    )

    diagnostic_columns <- c(
      "selected calibration level",
      "selected calibration level for x_1 = 0",
      "selected calibration level for x_1 = 1",
      "estimated calibration risk",
      "estimated calibration risk for x_1 = 0",
      "estimated calibration risk for x_1 = 1",
      "absolute calibration-to-test risk error",
      "absolute calibration-to-test risk error for x_1 = 0",
      "absolute calibration-to-test risk error for x_1 = 1",
      "fraction of censoring probabilities truncated at eta"
    )
    metrics[diagnostic_columns] <- NA_real_

    if (!is.null(calibration_diagnostics)) {
      for (row_idx in seq_len(nrow(metrics))) {
        method_name <- colnames(output_df)[row_idx]
        method_diag <- calibration_diagnostics[
          calibration_diagnostics$method == method_name,
          ,
          drop = FALSE
        ]

        if (nrow(method_diag) == 0L) next

        diag0 <- method_diag[method_diag$subgroup == 0, , drop = FALSE]
        diag1 <- method_diag[method_diag$subgroup == 1, , drop = FALSE]
        valid_weight <- is.finite(method_diag$n_calib) & method_diag$n_calib > 0

        if (any(valid_weight)) {
          metrics[["selected calibration level"]][row_idx] <- weighted.mean(
            method_diag$selected_calibration_level[valid_weight],
            method_diag$n_calib[valid_weight],
            na.rm = TRUE
          )
          metrics[["estimated calibration risk"]][row_idx] <- weighted.mean(
            method_diag$estimated_calibration_risk[valid_weight],
            method_diag$n_calib[valid_weight],
            na.rm = TRUE
          )
          metrics[["fraction of censoring probabilities truncated at eta"]][row_idx] <- weighted.mean(
            method_diag$fraction_censoring_probabilities_truncated[valid_weight],
            method_diag$n_calib[valid_weight],
            na.rm = TRUE
          )
        }

        if (nrow(diag0) == 1L) {
          metrics[["selected calibration level for x_1 = 0"]][row_idx] <- diag0$selected_calibration_level
          metrics[["estimated calibration risk for x_1 = 0"]][row_idx] <- diag0$estimated_calibration_risk
        }
        if (nrow(diag1) == 1L) {
          metrics[["selected calibration level for x_1 = 1"]][row_idx] <- diag1$selected_calibration_level
          metrics[["estimated calibration risk for x_1 = 1"]][row_idx] <- diag1$estimated_calibration_risk
        }

        metrics[["absolute calibration-to-test risk error"]][row_idx] <- abs(
          metrics[["estimated calibration risk"]][row_idx] -
            metrics[["true test-set miscoverage"]][row_idx]
        )
        metrics[["absolute calibration-to-test risk error for x_1 = 0"]][row_idx] <- abs(
          metrics[["estimated calibration risk for x_1 = 0"]][row_idx] -
            (1 - cov_grp0[row_idx])
        )
        metrics[["absolute calibration-to-test risk error for x_1 = 1"]][row_idx] <- abs(
          metrics[["estimated calibration risk for x_1 = 1"]][row_idx] -
            (1 - cov_grp1[row_idx])
        )
      }
    }

    metrics
  }

  df_cams <- compute_metrics(
    cams_res$output,
    times_vec = time_cams,
    suffix_label = NULL,
    calibration_diagnostics = cams_res$diagnostics
  )

  # local_df_cams <- compute_metrics(
  #   local_cams_res$output,
  #   times_vec = local_time_cams,
  #   suffix_label = NULL
  # )
  
  ########################################
  ## Compute & Bind Final Results
  ########################################
  if (!only_cams) {
    df_joint <- compute_metrics(
      res_joint$output,
      times_vec = res_joint$times,
      suffix_label = "(Joint)"
    )
    df_subgroup <- compute_metrics(
      output_subgroup,
      times_vec = times_subgroup,
      suffix_label = "(Subgroup)"
    )
    # Append CAMS to the final CSV output
    simu_out <- rbind(df_cams, df_joint, df_subgroup)
  } else {
    # simu_out <- rbind(df_cams, local_df_cams)
    simu_out <- df_cams
  }

  rownames(simu_out) <- NULL
  return(simu_out)
}

save_local_cams_info_allseeds <- function(cams_res,
                                          setting,
                                          seed) {

  diag_dir <- file.path("../local_cams_diagnostics")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

  lambda_file <- file.path(diag_dir, "selected_lambdas_all.csv")
  risk_file <- file.path(diag_dir, "selected_risks_all.csv")
  diagnostics_file <- file.path(diag_dir, "cell_diagnostics_all.csv")
  summary_file <- file.path(diag_dir, "local_cams_summary_all.csv")

  write_replace_seed <- function(new_df, file) {
    if (nrow(new_df) == 0) {
      return(invisible(NULL))
    }

    if (file.exists(file)) {
      old_df <- read.csv(file, stringsAsFactors = FALSE)

      # Remove old rows for this same setting and seed, so reruns do not duplicate
      old_df <- old_df[!(old_df$setting == setting & old_df$seed == seed), , drop = FALSE]

      common_cols <- union(names(old_df), names(new_df))

      for (cc in setdiff(common_cols, names(old_df))) old_df[[cc]] <- NA
      for (cc in setdiff(common_cols, names(new_df))) new_df[[cc]] <- NA

      old_df <- old_df[, common_cols, drop = FALSE]
      new_df <- new_df[, common_cols, drop = FALSE]

      out_df <- rbind(old_df, new_df)
    } else {
      out_df <- new_df
    }

    write.csv(out_df, file, row.names = FALSE)
  }

  # --------------------------------------------------
  # 1. Selected lambdas, long format
  # --------------------------------------------------
  lambda_rows <- list()

  if (!is.null(cams_res$selected_lambdas)) {
    for (group_name in names(cams_res$selected_lambdas)) {

      group_obj <- cams_res$selected_lambdas[[group_name]]

      if (is.null(group_obj)) next

      for (method_name in names(group_obj)) {

        lambda_vec <- group_obj[[method_name]]

        if (is.null(lambda_vec)) {
          lambda_rows[[length(lambda_rows) + 1]] <- data.frame(
            setting = setting,
            seed = seed,
            group = group_name,
            method = method_name,
            cell = NA_character_,
            lambda = NA_real_
          )
        } else {
          lambda_rows[[length(lambda_rows) + 1]] <- data.frame(
            setting = setting,
            seed = seed,
            group = group_name,
            method = method_name,
            cell = names(lambda_vec),
            lambda = as.numeric(lambda_vec),
            row.names = NULL
          )
        }
      }
    }
  }

  lambda_df <- if (length(lambda_rows) > 0) {
    do.call(rbind, lambda_rows)
  } else {
    data.frame()
  }

  write_replace_seed(lambda_df, lambda_file)

  # --------------------------------------------------
  # 2. Selected risks, long format
  # --------------------------------------------------
  risk_rows <- list()

  if (!is.null(cams_res$selected_risks)) {
    for (group_name in names(cams_res$selected_risks)) {

      group_obj <- cams_res$selected_risks[[group_name]]

      if (is.null(group_obj)) next

      for (method_name in names(group_obj)) {

        risk_obj <- group_obj[[method_name]]

        if (is.null(risk_obj)) {
          risk_rows[[length(risk_rows) + 1]] <- data.frame(
            setting = setting,
            seed = seed,
            group = group_name,
            method = method_name,
            risk_name = NA_character_,
            risk_value = NA_real_
          )

        } else {
          risk_vec <- unlist(risk_obj)

          risk_rows[[length(risk_rows) + 1]] <- data.frame(
            setting = setting,
            seed = seed,
            group = group_name,
            method = method_name,
            risk_name = names(risk_vec),
            risk_value = as.numeric(risk_vec),
            row.names = NULL
          )
        }
      }
    }
  }

  risk_df <- if (length(risk_rows) > 0) {
    do.call(rbind, risk_rows)
  } else {
    data.frame()
  }

  write_replace_seed(risk_df, risk_file)

  # --------------------------------------------------
  # 3. Cell diagnostics
  # --------------------------------------------------
  diag_rows <- list()

  if (!is.null(cams_res$diagnostics)) {
    for (group_name in names(cams_res$diagnostics)) {

      diag_df <- cams_res$diagnostics[[group_name]]

      if (is.null(diag_df) || nrow(diag_df) == 0) next

      diag_df$setting <- setting
      diag_df$seed <- seed
      diag_df$group <- group_name

      first_cols <- c("setting", "seed", "group")
      diag_df <- diag_df[, c(first_cols, setdiff(names(diag_df), first_cols)), drop = FALSE]

      diag_rows[[length(diag_rows) + 1]] <- diag_df
    }
  }

  diagnostics_df <- if (length(diag_rows) > 0) {
    do.call(rbind, diag_rows)
  } else {
    data.frame()
  }

  write_replace_seed(diagnostics_df, diagnostics_file)

  # --------------------------------------------------
  # 4. Summary file, one row per seed/group/method
  # --------------------------------------------------
  summary_rows <- list()

  if (!is.null(cams_res$selected_lambdas)) {
    for (group_name in names(cams_res$selected_lambdas)) {

      lambda_group <- cams_res$selected_lambdas[[group_name]]
      risk_group <- cams_res$selected_risks[[group_name]]

      if (is.null(lambda_group)) next

      for (method_name in names(lambda_group)) {

        lambda_vec <- lambda_group[[method_name]]

        lambda_x2_le <- NA_real_
        lambda_x2_gt <- NA_real_

        if (!is.null(lambda_vec)) {
          if ("X2<=0" %in% names(lambda_vec)) {
            lambda_x2_le <- as.numeric(lambda_vec["X2<=0"])
          }
          if ("X2>0" %in% names(lambda_vec)) {
            lambda_x2_gt <- as.numeric(lambda_vec["X2>0"])
          }
        }

        group_risk <- NA_real_
        cell_risk_le <- NA_real_
        cell_risk_gt <- NA_real_

        if (!is.null(risk_group) && method_name %in% names(risk_group)) {
          risk_obj <- risk_group[[method_name]]

          if (!is.null(risk_obj)) {
            risk_vec <- unlist(risk_obj)

            if ("group_risk" %in% names(risk_vec)) {
              group_risk <- as.numeric(risk_vec["group_risk"])
            }

            if ("cell_risk_X2<=0" %in% names(risk_vec)) {
              cell_risk_le <- as.numeric(risk_vec["cell_risk_X2<=0"])
            }

            if ("cell_risk_X2>0" %in% names(risk_vec)) {
              cell_risk_gt <- as.numeric(risk_vec["cell_risk_X2>0"])
            }
          }
        }

        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          setting = setting,
          seed = seed,
          group = group_name,
          method = method_name,
          lambda_X2_le = lambda_x2_le,
          lambda_X2_gt = lambda_x2_gt,
          group_risk = group_risk,
          cell_risk_X2_le = cell_risk_le,
          cell_risk_X2_gt = cell_risk_gt
        )
      }
    }
  }

  summary_df <- if (length(summary_rows) > 0) {
    do.call(rbind, summary_rows)
  } else {
    data.frame()
  }

  write_replace_seed(summary_df, summary_file)

  invisible(list(
    lambda_file = lambda_file,
    risk_file = risk_file,
    diagnostics_file = diagnostics_file,
    summary_file = summary_file
  ))
}
