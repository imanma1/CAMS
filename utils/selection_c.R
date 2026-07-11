selection_c <- function(data, p, n, xnames,
                        c_ref, weight_ref,
                        model = "cox",
                        type = "quantile",
                        dist = "weibull",
                        mdl0 = NULL,
                        alpha,
                        use_oracle_sc = FALSE) {

  library(parallel)

  if (length(c_ref) == 0) {
    return(list(
      c_opt = NA_real_,
      c_ref = c_ref,
      bnd_ref = numeric(0),
      selection_failed = TRUE
    ))
  }

  if (.Platform$OS.type == "windows") {
    n_threads <- 1L
  } else {
    slurm_cores <- suppressWarnings(
      as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
    )

    if (is.na(slurm_cores) || slurm_cores < 1L) {
      n_threads <- parallel::detectCores()
    } else {
      n_threads <- slurm_cores
    }

    n_threads <- min(n_threads, length(c_ref))
  }

  if (.Platform$OS.type != "windows") {
    library(RhpcBLASctl)

    original_blas_threads <- blas_get_num_procs()

    blas_set_num_threads(1)
    omp_set_num_threads(1)

    on.exit(
      {
        blas_set_num_threads(original_blas_threads)
        omp_set_num_threads(original_blas_threads)
      },
      add = TRUE
    )
  }

  bnd_ref <- unlist(
    mclapply(
      seq_along(c_ref),
      function(i) {

        current_weight <- if (is.null(weight_ref)) {
          NULL
        } else {
          weight_ref[, i]
        }

        tryCatch(
          evaluate_length(
            c = c_ref[i],
            alpha = alpha,
            n = n,
            p = p,
            model = model,
            data = data,
            weight = current_weight,
            xnames = xnames,
            type = type,
            dist = dist,
            mdl0 = mdl0,
            use_oracle_sc = use_oracle_sc
          ),
          error = function(e) {
            -Inf
          }
        )
      },
      mc.cores = n_threads
    )
  )

  bnd_ref <- as.numeric(bnd_ref)

  # Ensure one result exists for every cutoff candidate.
  if (length(bnd_ref) != length(c_ref)) {
    corrected_bnd_ref <- rep(-Inf, length(c_ref))

    n_copy <- min(
      length(bnd_ref),
      length(c_ref)
    )

    if (n_copy > 0) {
      corrected_bnd_ref[seq_len(n_copy)] <-
        bnd_ref[seq_len(n_copy)]
    }

    bnd_ref <- corrected_bnd_ref
  }

  valid_idx <- which(
    is.finite(c_ref) &
      is.finite(bnd_ref)
  )

  cat(
    sprintf(
      "Valid cutoff candidates: %d/%d\n",
      length(valid_idx),
      length(c_ref)
    )
  )

  if (length(valid_idx) == 0) {
    warning(
      paste0(
        "selection_c(): all cutoff candidates failed. ",
        "The subgroup fitting sample may be too small ",
        "relative to the number of predictors."
      )
    )

    return(list(
      c_opt = NA_real_,
      c_ref = c_ref,
      bnd_ref = bnd_ref,
      selection_failed = TRUE
    ))
  }

  best_idx <- valid_idx[
    which.max(bnd_ref[valid_idx])
  ]

  c_opt <- c_ref[best_idx]

  return(list(
    c_opt = c_opt,
    c_ref = c_ref,
    bnd_ref = bnd_ref,
    selection_failed = FALSE
  ))
}


evaluate_length <- function(c, alpha, n, p,
                            model, data, weight, xnames,
                            type = "quantile",
                            dist = "weibull",
                            seed = 2020,
                            mdl0 = NULL,
                            use_oracle_sc = FALSE) {

  if (
    length(c) != 1L ||
    is.na(c) ||
    !is.finite(c)
  ) {
    return(-Inf)
  }

  set.seed(seed)

  if (n < 4L) {
    return(-Inf)
  }

  I_fit <- sample(
    seq_len(n),
    floor(n / 2),
    replace = FALSE
  )

  remaining_idx <- setdiff(
    seq_len(n),
    I_fit
  )

  I_calib <- sample(
    remaining_idx,
    floor(n / 4),
    replace = FALSE
  )

  I_test <- setdiff(
    remaining_idx,
    I_calib
  )

  data_fit <- data[I_fit, , drop = FALSE]
  data_calib <- data[I_calib, , drop = FALSE]
  data_test <- data[I_test, , drop = FALSE]

  # survreg cannot meaningfully estimate more regression
  # coefficients than the number of internal fitting rows.
  if (nrow(data_fit) <= length(xnames) + 1L) {
    return(-Inf)
  }

  if (
    nrow(data_calib) == 0 ||
    nrow(data_test) == 0
  ) {
    return(-Inf)
  }

  if (is.null(weight)) {

    censoring_probabilities <- tryCatch(
      {
        if (
          use_oracle_sc ||
          inherits(mdl0, "oracle_sc")
        ) {

          list(
            pr_calib = sc_prob(
              mdl0 = mdl0,
              data = data_calib,
              xnames = xnames,
              t = c
            ),
            pr_new = sc_prob(
              mdl0 = mdl0,
              data = data_test,
              xnames = xnames,
              t = c
            )
          )

        } else {

          cens_prob(
            mdl = mdl0,
            calib = data_calib,
            test = data_test,
            method = "gpr",
            xnames = xnames,
            c = c
          )
        }
      },
      error = function(e) {
        NULL
      }
    )

    if (is.null(censoring_probabilities)) {
      return(-Inf)
    }

    pr_calib <- censoring_probabilities$pr_calib
    pr_new <- censoring_probabilities$pr_new

    if (
      length(pr_calib) != nrow(data_calib) ||
      length(pr_new) != nrow(data_test) ||
      any(!is.finite(pr_calib)) ||
      any(!is.finite(pr_new)) ||
      any(pr_calib <= 0) ||
      any(pr_new <= 0)
    ) {
      return(-Inf)
    }

    weight_calib <- 1 / pr_calib
    weight_new <- 1 / pr_new

  } else {

    if (length(weight) != nrow(data)) {
      return(-Inf)
    }

    weight_calib <- weight[I_calib]
    weight_new <- weight[I_test]
  }

  if (
    any(!is.finite(weight_calib)) ||
    any(!is.finite(weight_new))
  ) {
    return(-Inf)
  }

  x <- data_test[, xnames, drop = FALSE]

  bnd <- tryCatch(
    {
      cox0_based(
        x = x,
        p = p,
        len_x = nrow(x),
        xnames = xnames,
        c = c,
        alpha = alpha,
        data_fit = data_fit,
        data_calib = data_calib,
        type = type,
        dist = dist,
        weight_calib = weight_calib,
        weight_new = weight_new
      )
    },
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(bnd) ||
    length(bnd) != nrow(x) ||
    any(!is.finite(bnd))
  ) {
    return(-Inf)
  }

  mean_bnd <- mean(bnd)

  if (!is.finite(mean_bnd)) {
    return(-Inf)
  }

  return(mean_bnd)
}