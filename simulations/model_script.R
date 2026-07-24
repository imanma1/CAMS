draw_weibull_aft_time <- function(
    mu_t,
    sigma_t = 0.45
) {

  n <- length(mu_t)

  if (length(sigma_t) == 1L) {
    sigma_t <- rep(sigma_t, n)
  }

  if (length(sigma_t) != n) {
    stop("sigma_t must have length 1 or the same length as mu_t.")
  }

  u <- pmin(
    pmax(runif(n), 1e-12),
    1 - 1e-12
  )

  eps_t <- log(-log(u))

  exp(
    mu_t +
      sigma_t * eps_t
  )
}


shared_weibull_mu_20 <- function(x) {

  x <- as.data.frame(x)

  required_names <- paste0("X", 1:20)

  missing_names <- setdiff(
    required_names,
    colnames(x)
  )

  if (length(missing_names) > 0L) {
    stop(
      sprintf(
        "The shared 20-dimensional event model is missing: %s",
        paste(missing_names, collapse = ", ")
      )
    )
  }

  # Coefficients 0.05 * (-1)^j for j = 5, ..., 20.
  beta_dense <- 0.05 * (-1)^(5:20)

  dense_score <- as.numeric(
    as.matrix(
      x[, paste0("X", 5:20), drop = FALSE]
    ) %*% beta_dense
  )

  2.40 +
    0.50 * x$X1 +
    0.70 * x$X2 -
    0.50 * x$X3 +
    0.30 * x$X4 +
    dense_score
}


basis_shared_weibull_mu <- function(x) {

  x <- as.data.frame(x)

  2.20 +
    0.60 * x$X2 -
    0.40 * x$X3 +
    0.25 * x$X4 +
    0.30 * x$X5 +
    0.45 * x$X1
}

draw_min_extreme_value_time <- function(
    mu_t,
    sigma_t
) {

  n <- length(mu_t)

  if (length(sigma_t) == 1L) {
    sigma_t <- rep(sigma_t, n)
  }

  if (length(sigma_t) != n) {
    stop(
      paste0(
        "sigma_t must have length one or ",
        "the same length as mu_t."
      )
    )
  }

  u <- pmin(
    pmax(
      runif(n),
      1e-12
    ),
    1 - 1e-12
  )

  eps_t <- log(
    -log(u)
  )

  exp(
    mu_t +
      sigma_t * eps_t
  )
}


shared_hd_event_location <- function(x) {

  x <- as.data.frame(x)

  required_names <- paste0(
    "X",
    1:75
  )

  missing_names <- setdiff(
    required_names,
    colnames(x)
  )

  if (length(missing_names) > 0L) {
    stop(
      sprintf(
        "The shared HD event model is missing: %s",
        paste(
          missing_names,
          collapse = ", "
        )
      )
    )
  }

  beta_dense <- 0.03 * rep(
    c(1, -1),
    length.out = 75 - 4
  )

  dense_score <- as.numeric(
    as.matrix(
      x[
        ,
        paste0("X", 5:75),
        drop = FALSE
      ]
    ) %*% beta_dense
  )

  2.60 +
    0.80 * x$X1 +
    0.60 * x$X2 -
    0.50 * x$X3 +
    0.30 * x$X4 +
    dense_score
}


mild_intersection_censoring_location <- function(x) {

  x <- as.data.frame(x)

  3.25 -
    0.15 * x$X1 -
    0.10 * as.numeric(x$X2 > 0) +
    0.15 * x$X3
}

model_generating_fun <- function(n_train, n_calib, n_test,
                                 setting, xmin, xmax,
                                 bernoulli_prob = 0.1,
                                 homoscedastic_event = FALSE) {

  # =====================================================================
  # 1. DEFINE SETTING FORMULAS AND DIMENSIONS
  # =====================================================================
  if (setting == "homo_cens") {
    # Setting (i): Homogeneous censoring
    p <- 2
    gen_t <- function(x) exp(2 + 0.5 * x[,2] - 0.5 * x[,1] + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) exp(3.0 + 0.5 * rnorm(nrow(x)))

  } else if (setting == "cov_cens") {
    # Setting (ii): Censoring depending on ordinary covariates
    p <- 2
    gen_t <- function(x) exp(2 + 0.5 * x[,2] - 0.5 * x[,1] + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) exp(2.5 + 1.2 * x[,2] + 0.5 * rnorm(nrow(x)))

  } else if (setting == "prot_cens") {
    # Setting (iii): Censoring depending on Protected Group
    p <- 2
    gen_t <- function(x) exp(2 + 0.5 * x[,2] - 0.5 * x[,1] + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) exp(3.0 - 1.5 * x[,1] + 0.5 * rnorm(nrow(x)))

  } else if (setting == "heavy_prot_cens") {
    # Setting (iv): Heavy censoring in one protected group
    p <- 2
    gen_t <- function(x) exp(2 + 0.5 * x[,2] + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) exp(3.5 - 2.8 * x[,1] + 0.5 * rnorm(nrow(x)))

  } else if (setting == "heavy_inter_cens") {
    # Setting (v): Heavy censoring in an intersectional subgroup
    p <- 2
    gen_t <- function(x) exp(2 + 0.5 * x[,2] + 0.5 * rnorm(nrow(x)))
    # R automatically coerces the boolean (x[,1]==1 & x[,2]>0) into 1s and 0s
    gen_c <- function(x) exp(3.5 - 1.8 * (x[,1] == 1 & x[,2] > 0) + 0.5 * rnorm(nrow(x)))

  } else if (setting == "surv_misspec") {
    # Setting (vi): Survival-model misspecification
    p <- 3
    gen_t <- function(x) exp(2 + 0.5 * (x[,2]^2) + sin(x[,3]) * x[,1] + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) exp(2.5 + 0.5 * x[,2] + 0.5 * rnorm(nrow(x)))

  } else if (setting == "cens_misspec") {
    # Setting (vii): Censoring-model misspecification (10 Variables)
    p <- 10
    gen_t <- function(x) exp(2 + 0.5 * x[,2] - 0.5 * x[,1] + 0.1 * rowSums(x[, 4:10, drop=FALSE]) + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) exp(2.5 + (x[,2]^2) * x[,1] + cos(x[,3]) + 0.2 * rowSums(x[, 4:10, drop=FALSE]^2) + 0.5 * rnorm(nrow(x)))

  } else if (setting == "simul_misspec") {
    # Setting (viii): Simultaneous moderate misspecification (10 Variables)
    p <- 10
    # X_i * X_{i-1} is calculated via vectorized element-wise matrix multiplication
    gen_t <- function(x) exp(2 + 0.5 * (x[,2]^2) + 0.5 * x[,2] * x[,1] + 
                               0.1 * rowSums(x[, 3:10, drop=FALSE] * x[, 2:9, drop=FALSE]) + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) exp(2.5 + (x[,3]^2) - 1.5 * x[,1] + 
                               0.1 * rowSums(abs(x[, 4:10, drop=FALSE])) + 0.5 * rnorm(nrow(x)))
                               
  } else if (setting == "complex_surv") {
    # Setting: Complex nonlinear survival with homogeneous censoring (10 Variables)
    p <- 10
    # T = exp(2 - 1.5*X1 + 0.5*X2^2 + 0.3*X3*X4 + 0.1*sum(X5:X10) + 0.5*eps)
    gen_t <- function(x) {
      exp(2 - 1.5 * x[,1] + 0.5 * (x[,2]^2) + 0.3 * x[,3] * x[,4] + 
          0.1 * rowSums(x[, 5:10, drop=FALSE]) + 0.5 * rnorm(nrow(x)))
    }

    # C = exp(3.0 + 0.5*eps)
    gen_c <- function(x) {
      exp(3.0 + 0.5 * rnorm(nrow(x)))
    }
  } else if (setting == "var_shift_heavy_cens") {
    # Setting: Variance Shift + Heavy Censoring (10 Variables)
    p <- 10

    # T = exp(2 + 0.5*X2 - 1.0*X1 + (0.5 + 0.8*X1)*eps_T)
    gen_t <- function(x) {
      exp(2 + 0.5 * x[,2] - 1.0 * x[,1] + 
         (0.5 + 0.8 * x[,1]) * rnorm(nrow(x)))
    }

    # C = exp(3.0 + 0.5*X2 - 2.8*X1 + 0.5*eps_C)
    gen_c <- function(x) {
      exp(3.0 + 0.5 * x[,2] - 2.8 * x[,1] + 
          0.5 * rnorm(nrow(x)))
    }
  } else if (setting == "starve_hetero") {
    # Setting: Starvation + Heteroscedasticity (15 Variables)
    p <- 15

    # T = exp(2 - 1.5*X1 + sin(X2*X3) + 0.5*X4^2 + 0.2*sum(X5:X15) + (0.5 + 0.8*X1)*eps_T)
    gen_t <- function(x) {
      exp(2 - 1.5 * x[,1] + sin(x[,2] * x[,3]) + 0.5 * (x[,4]^2) + 
          0.2 * rowSums(x[, 5:15, drop=FALSE]) + 
          (0.5 + 0.8 * x[,1]) * rnorm(nrow(x)))
    }

    # C = exp(3.0 - 1.5*X1 + 0.5*eps_C)
    gen_c <- function(x) {
      exp(3.0 - 1.5 * x[,1] + 0.5 * rnorm(nrow(x)))
    }

  } else if (setting == "starve_hetero_high_dim") {
    # Setting: Starvation + Heteroscedasticity (75 Variables)
    p <- 75

    # T = exp(2 - 1.5*X1 + 0.1*sum(X2:X75) + (0.5 + 0.8*X1)*eps_T)
    gen_t <- function(x) {
      exp(2 - 1.5 * x[,1] + 
          0.1 * rowSums(x[, 2:75, drop=FALSE]) + 
          (0.5 + 0.8 * x[,1]) * rnorm(nrow(x)))
    }

    # C = exp(3.0 - 1.5*X1 + 0.5*eps_C)
    gen_c <- function(x) {
      exp(3.0 - 1.5 * x[,1] + 0.5 * rnorm(nrow(x)))
    }

  } else if (setting == "high_survival_heavy_cens") {

    # Long-survival protected group, but that group is more heavily censored
    p <- 3

    # X1 = 1 has genuinely longer survival
    gen_t <- function(x) {
      exp(
        2.0 +
        0.60 * x[, 2] -
        0.30 * x[, 3] +
        0.70 * x[, 1] +
        0.50 * rnorm(nrow(x))
      )
    }

    # X1 = 1 is censored earlier.
    # High X2 also means higher T but earlier C.
    gen_c <- function(x) {
      exp(
        3.40 -
        0.50 * x[, 2] +
        0.20 * x[, 3] -
        0.80 * x[, 1] +
        0.50 * rnorm(nrow(x))
      )
    }
  } else if (setting == "anti_aligned_cens") {

    p <- 5

    gen_t <- function(x) {
      sigma_t <- 0.45 + 0.20 * x[, 1]

      exp(
        2.0 +
        0.70 * x[, 2] -
        0.50 * x[, 3] +
        0.25 * x[, 4] -
        0.40 * x[, 1] +
        sigma_t * rnorm(nrow(x))
      )
    }

    gen_c <- function(x) {
      exp(
        2.80 -
        0.70 * x[, 2] +
        0.50 * x[, 3] -
        0.80 * x[, 1] +
        0.50 * rnorm(nrow(x))
      )
    }
  } else if (setting == "weibull_aft_anti_cens") {

    p <- 5

    gen_t <- function(x) {

      # Minimum extreme-value error used by a Weibull AFT model
      eps_t <- log(-log(runif(nrow(x))))

      exp(
        2.0 +
        0.70 * x[, 2] -
        0.50 * x[, 3] +
        0.30 * x[, 4] +
        0.60 * x[, 1] +
        0.45 * eps_t
      )
    }

    gen_c <- function(x) {
      exp(
        3.30 -
        0.60 * x[, 2] +
        0.40 * x[, 3] -
        0.80 * x[, 1] +
        0.50 * rnorm(nrow(x))
      )
    } 
  } else if (setting == "moderate_inter_cens") {

    p <- 5

    gen_t <- function(x) {
      exp(
        2.0 +
        0.65 * x[, 2] -
        0.45 * x[, 3] +
        0.25 * x[, 4] -
        0.30 * x[, 1] +
        0.50 * rnorm(nrow(x))
      )
    }

    hard_region <- function(x) {
      (x[, 2] > 0) & (x[, 3] < 0)
    }

    gen_c <- function(x) {
      exp(
        3.00 -
        0.40 * hard_region(x) -
        0.60 * x[, 1] +
        0.20 * x[, 4] +
        0.50 * rnorm(nrow(x))
      )
    }
  } else if (setting == "cams_pool_calib_hd") {

    # Rare minority + dense shared high-dimensional signal
    # + group-specific survival scale
    p <- 75

    # All 74 ordinary covariates contribute to survival.
    # The signal is shared by X1 = 0 and X1 = 1.
    beta <- 0.05 * rep(c(1, -1), length.out = p - 1)

    dense_score <- function(x) {
      as.numeric(
        as.matrix(x[, 2:p, drop = FALSE]) %*% beta
      )
    }

    gen_t <- function(x) {

      score <- dense_score(x)

      # Standard extreme-value noise, matching survreg(dist = "weibull")
      u <- runif(nrow(x))
      u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
      eps_t <- log(-log(u))

      # Majority has a narrow survival distribution.
      # Minority has a substantially wider survival distribution.
      sigma_t <- ifelse(x[, 1] == 1, 0.95, 0.35)

      # The positive minority intercept keeps the minority's lower tail
      # scientifically informative despite its greater variance.
      mu_t <- 2.2 +
        1.35 * x[, 1] +
        score

      exp(mu_t + sigma_t * eps_t)
    }

    gen_c <- function(x) {

      score <- dense_score(x)

      # Moderate rather than nearly complete censoring.
      # Censoring is somewhat heavier in X1 = 1 and occurs earlier
      # in regions where true survival tends to be longer.
      exp(
        3.15 -
        0.15 * x[, 1] -
        0.10 * score +
        0.50 * rnorm(nrow(x))
      )
    }
  } else if (setting == "cams_pool_calib_hd_mild") {

    p <- 75

    beta <- 0.05 * rep(c(1, -1), length.out = p - 1)

    dense_score <- function(x) {
      as.numeric(
        as.matrix(x[, 2:p, drop = FALSE]) %*% beta
      )
    }

    gen_t <- function(x) {

      score <- dense_score(x)

      u <- pmin(pmax(runif(nrow(x)), 1e-12), 1 - 1e-12)
      eps_t <- log(-log(u))

      sigma_t <- ifelse(x[, 1] == 1, 0.85, 0.35)

      mu_t <- 2.2 +
        1.13 * x[, 1] +
        score

      exp(mu_t + sigma_t * eps_t)
    }

    gen_c <- function(x) {

      score <- dense_score(x)

      exp(
        3.20 -
        0.10 * x[, 1] -
        0.08 * score +
        0.50 * rnorm(nrow(x))
      )
    }
  } else if (setting == "cams_pool_calib_hd_strong") {

    p <- 75

    beta <- 0.05 * rep(c(1, -1), length.out = p - 1)

    dense_score <- function(x) {
      as.numeric(
        as.matrix(x[, 2:p, drop = FALSE]) %*% beta
      )
    }

    gen_t <- function(x) {

      score <- dense_score(x)

      u <- pmin(pmax(runif(nrow(x)), 1e-12), 1 - 1e-12)
      eps_t <- log(-log(u))

      sigma_t <- ifelse(x[, 1] == 1, 1.05, 0.35)

      mu_t <- 2.2 +
        1.58 * x[, 1] +
        score

      exp(mu_t + sigma_t * eps_t)
    }

    gen_c <- function(x) {

      score <- dense_score(x)

      exp(
        3.20 -
        0.10 * x[, 1] -
        0.10 * score +
        0.50 * rnorm(nrow(x))
      )
    }
  } else if (setting == "cams_vs_vanilla_lower_tail_hd_mild") {

    # Rare subgroup + high-dimensional shared survival pattern
    p <- 75

    beta_dense <- 0.03 * rep(
      c(1, -1),
      length.out = p - 4
    )

    dense_score <- function(x) {
      as.numeric(
        as.matrix(x[, 5:p, drop = FALSE]) %*% beta_dense
      )
    }

    mu_t_fun <- function(x) {
      2.8 +
        0.40 * x[, 1] +
        0.60 * x[, 2] -
        0.50 * x[, 3] +
        0.30 * x[, 4] +
        dense_score(x)
    }

    sigma_t_fun <- function(x) {
      if (homoscedastic_event) {
        rep(0.35, nrow(x))
      } else {
        0.28 + 0.17 * x[, 1]
      }
    }

    gen_t <- function(x) {

      # Extreme-value error:
      # log(T) follows a Weibull AFT model compatible with survreg.
      u <- pmin(
        pmax(runif(nrow(x)), 1e-12),
        1 - 1e-12
      )

      eps_t <- log(-log(u))

      exp(
        mu_t_fun(x) +
          sigma_t_fun(x) * eps_t
      )
    }

    gen_c <- function(x) {

      mu_t <- mu_t_fun(x)

      # Early censoring occurs for 15% of the majority
      # and 20% of the minority.
      prob_early <- 0.10 + 0.05 * x[, 1]

      early <- rbinom(
        nrow(x),
        size = 1,
        prob = prob_early
      )

      # Early component lies near/below the target lower tail.
      mu_early <- mu_t - 1.00

      # Late component usually lies well beyond the event time.
      mu_late <- mu_t + 1.20

      log_c <- ifelse(
        early == 1,
        mu_early + 0.25 * rnorm(nrow(x)),
        mu_late + 0.40 * rnorm(nrow(x))
      )

      exp(log_c)
    }
  } else if (setting == "cams_vs_vanilla_lower_tail_hd_main") {

    # Rare subgroup + high-dimensional shared survival pattern
    p <- 75

    beta_dense <- 0.03 * rep(
      c(1, -1),
      length.out = p - 4
    )

    dense_score <- function(x) {
      as.numeric(
        as.matrix(x[, 5:p, drop = FALSE]) %*% beta_dense
      )
    }

    mu_t_fun <- function(x) {
      2.8 +
        0.40 * x[, 1] +
        0.60 * x[, 2] -
        0.50 * x[, 3] +
        0.30 * x[, 4] +
        dense_score(x)
    }

    sigma_t_fun <- function(x) {
      if (homoscedastic_event) {
        rep(0.35, nrow(x))
      } else {
        0.28 + 0.17 * x[, 1]
      }
    }

    gen_t <- function(x) {

      # Extreme-value error:
      # log(T) follows a Weibull AFT model compatible with survreg.
      u <- pmin(
        pmax(runif(nrow(x)), 1e-12),
        1 - 1e-12
      )

      eps_t <- log(-log(u))

      exp(
        mu_t_fun(x) +
          sigma_t_fun(x) * eps_t
      )
    }

    gen_c <- function(x) {

      mu_t <- mu_t_fun(x)

      # Early censoring occurs for 15% of the majority
      # and 20% of the minority.
      prob_early <- 0.15 + 0.05 * x[, 1]

      early <- rbinom(
        nrow(x),
        size = 1,
        prob = prob_early
      )

      # Early component lies near/below the target lower tail.
      mu_early <- mu_t - 1.30

      # Late component usually lies well beyond the event time.
      mu_late <- mu_t + 1.20

      log_c <- ifelse(
        early == 1,
        mu_early + 0.25 * rnorm(nrow(x)),
        mu_late + 0.40 * rnorm(nrow(x))
      )

      exp(log_c)
    }
  } else if (setting == "cams_vs_vanilla_lower_tail_hd_strong") {

    # Rare subgroup + high-dimensional shared survival pattern
    p <- 75

    beta_dense <- 0.03 * rep(
      c(1, -1),
      length.out = p - 4
    )

    dense_score <- function(x) {
      as.numeric(
        as.matrix(x[, 5:p, drop = FALSE]) %*% beta_dense
      )
    }

    mu_t_fun <- function(x) {
      2.8 +
        0.40 * x[, 1] +
        0.60 * x[, 2] -
        0.50 * x[, 3] +
        0.30 * x[, 4] +
        dense_score(x)
    }

    sigma_t_fun <- function(x) {
      if (homoscedastic_event) {
        rep(0.35, nrow(x))
      } else {
        0.28 + 0.17 * x[, 1]
      }
    }

    gen_t <- function(x) {

      # Extreme-value error:
      # log(T) follows a Weibull AFT model compatible with survreg.
      u <- pmin(
        pmax(runif(nrow(x)), 1e-12),
        1 - 1e-12
      )

      eps_t <- log(-log(u))

      exp(
        mu_t_fun(x) +
          sigma_t_fun(x) * eps_t
      )
    }

    gen_c <- function(x) {

      mu_t <- mu_t_fun(x)

      # Early censoring occurs for 15% of the majority
      # and 20% of the minority.
      prob_early <- 0.20 + 0.05 * x[, 1]

      early <- rbinom(
        nrow(x),
        size = 1,
        prob = prob_early
      )

      # Early component lies near/below the target lower tail.
      mu_early <- mu_t - 1.50

      # Late component usually lies well beyond the event time.
      mu_late <- mu_t + 1.20

      log_c <- ifelse(
        early == 1,
        mu_early + 0.25 * rnorm(nrow(x)),
        mu_late + 0.40 * rnorm(nrow(x))
      )

      exp(log_c)
    }
  } else if (
    setting %in% c(
      "local_null_constant_scale",
      "local_within_group_scale_x2",
      "local_interaction_scale_x2_x3",
      "local_unsupported_rare_pocket",
      "local_oblique_scale_x2_x3"
    )
  ) {

    # ============================================================
    # Local-CAMS testing settings
    # ============================================================
    p <- 75

    beta_dense <- 0.03 * rep(
      c(1, -1),
      length.out = p - 4
    )

    dense_score <- function(x) {
      as.numeric(
        as.matrix(
          x[, 5:p, drop = FALSE]
        ) %*% beta_dense
      )
    }

    # Same event-time location model in every setting
    mu_t_fun <- function(x) {
      2.8 +
        0.40 * x[, 1] +
        0.60 * x[, 2] -
        0.50 * x[, 3] +
        0.30 * x[, 4] +
        dense_score(x)
    }

    # Setting-specific event-time scale
    sigma_t_fun <- function(x) {

      if (setting == "local_null_constant_scale") {

        # Negative control:
        # no within-R lower-tail heterogeneity
        sigma_t <- rep(
          0.35,
          nrow(x)
        )

      } else if (setting == "local_within_group_scale_x2") {

        # Simple positive control:
        # useful split should be near X2 = 0
        sigma_t <- 0.25 +
          0.16 * as.numeric(x[, 2] > 0) +
          0.05 * x[, 1]

      } else if (setting == "local_interaction_scale_x2_x3") {

        # Harder positive control:
        # high-scale region requires an interaction
        sigma_t <- 0.25 +
          0.18 * as.numeric(
            x[, 2] > 0 &
              x[, 3] > 0
          ) +
          0.05 * x[, 1]

      } else if (setting == "local_unsupported_rare_pocket") {

        # Rare region should generally not be split because
        # it lacks sufficient audit/calibration support
        sigma_t <- 0.25 +
          0.22 * as.numeric(
            x[, 2] > 1.4 &
              x[, 3] > 1.4
          ) +
          0.05 * x[, 1]

      } else if (setting == "local_oblique_scale_x2_x3") {

        # Non-axis-aligned positive control.
        # Useful for testing the limitation of a shallow
        # axis-aligned regression tree.
        sigma_t <- 0.25 +
          0.18 * as.numeric(
            0.80 * x[, 2] -
              0.60 * x[, 3] > 0
          ) +
          0.05 * x[, 1]
      }

      as.numeric(sigma_t)
    }

    gen_t <- function(x) {

      u <- pmin(
        pmax(
          runif(nrow(x)),
          1e-12
        ),
        1 - 1e-12
      )

      eps_t <- log(
        -log(u)
      )

      exp(
        mu_t_fun(x) +
          sigma_t_fun(x) * eps_t
      )
    }

    # Same censoring DGP as
    # cams_vs_vanilla_lower_tail_hd_main
    gen_c <- function(x) {

      mu_t <- mu_t_fun(x)

      prob_early <- 0.15 +
        0.05 * x[, 1]

      early <- rbinom(
        nrow(x),
        size = 1,
        prob = prob_early
      )

      mu_early <- mu_t - 1.30
      mu_late <- mu_t + 1.20

      log_c <- ifelse(
        early == 1,
        mu_early +
          0.25 * rnorm(nrow(x)),
        mu_late +
          0.40 * rnorm(nrow(x))
      )

      exp(log_c)
    }
  } else if (
    setting == "rare_intersection_shared_weibull"
  ) {

    # ------------------------------------------------------------
    # Rare intersectional groups with a shared 20-dimensional
    # Weibull AFT event model.
    #
    # Intended configuration:
    # bernoulli_prob = 0.10
    # ------------------------------------------------------------

    p <- 20

    mu_t_fun <- function(x) {
      shared_weibull_mu_20(x)
    }

    gen_t <- function(x) {

      draw_weibull_aft_time(
        mu_t = mu_t_fun(x),
        sigma_t = 0.45
      )
    }

    gen_c <- function(x) {

      group_1_positive <- as.numeric(
        x[, 1] == 1 &
          x[, 2] > 0
      )

      group_1_nonpositive <- as.numeric(
        x[, 1] == 1 &
          x[, 2] <= 0
      )

      group_0_positive <- as.numeric(
        x[, 1] == 0 &
          x[, 2] > 0
      )

      mu_c <- 3.20 +
        0.20 * x[, 3] -
        0.60 * group_1_positive -
        0.30 * group_1_nonpositive -
        0.20 * group_0_positive

      exp(
        mu_c +
          0.50 * rnorm(nrow(x))
      )
    }
  } else if (
    setting == "basis_intersection_shared_weibull"
  ) {

    # ------------------------------------------------------------
    # The nonlinear structure is supplied to the fitted models as:
    #
    # X4 = X2^2
    # X5 = sin(X3)
    #
    # The pooled event model is therefore a correctly specified
    # linear Weibull AFT model in X1, ..., X5.
    # ------------------------------------------------------------

    p <- 5

    mu_t_fun <- function(x) {
      basis_shared_weibull_mu(x)
    }

    gen_t <- function(x) {

      draw_weibull_aft_time(
        mu_t = mu_t_fun(x),
        sigma_t = 0.45
      )
    }

    gen_c <- function(x) {

      group_1_positive <- as.numeric(
        x[, 1] == 1 &
          x[, 2] > 0
      )

      group_0_positive <- as.numeric(
        x[, 1] == 0 &
          x[, 2] > 0
      )

      mu_c <- 3.00 -
        0.65 * group_1_positive -
        0.35 * group_0_positive +
        0.25 * x[, 3]

      exp(
        mu_c +
          0.50 * rnorm(nrow(x))
      )
    }
  } else if (
    setting == "mixture_intersection_shared_weibull"
  ) {

    # ------------------------------------------------------------
    # Shared correctly specified Weibull AFT event model, with
    # censoring generated from an early/late lognormal mixture.
    #
    # The mixture depends on intersectional group membership.
    # ------------------------------------------------------------

    p <- 20

    mu_t_fun <- function(x) {
      shared_weibull_mu_20(x)
    }

    gen_t <- function(x) {

      draw_weibull_aft_time(
        mu_t = mu_t_fun(x),
        sigma_t = 0.45
      )
    }

    gen_c <- function(x) {

      n <- nrow(x)
      mu_t <- mu_t_fun(x)

      group_1_positive <- as.numeric(
        x[, 1] == 1 &
          x[, 2] > 0
      )

      group_0_positive <- as.numeric(
        x[, 1] == 0 &
          x[, 2] > 0
      )

      prob_early <- 0.08 +
        0.25 * group_1_positive +
        0.12 * group_0_positive

      prob_early <- pmin(
        pmax(prob_early, 0),
        1
      )

      early_component <- rbinom(
        n = n,
        size = 1,
        prob = prob_early
      )

      early_log_c <- (
        mu_t -
          0.35 +
          0.30 * rnorm(n)
      )

      late_log_c <- (
        mu_t +
          1.00 +
          0.40 * rnorm(n)
      )

      log_c <- ifelse(
        early_component == 1,
        early_log_c,
        late_log_c
      )

      exp(log_c)
    }
  } else if (
    setting == "intersection_scale_shift_hd"
  ) {

    # ============================================================
    # Shared high-dimensional event location, but four different
    # lower-tail scales across the intersectional R groups.
    #
    # R1: X1 = 0, X2 <= 0  -> sigma = 0.30
    # R2: X1 = 0, X2 >  0  -> sigma = 0.45
    # R3: X1 = 1, X2 <= 0  -> sigma = 0.55
    # R4: X1 = 1, X2 >  0  -> sigma = 0.85
    #
    # The censoring model is deliberately only mildly
    # heterogeneous so the main difficulty is event-tail
    # heterogeneity rather than positivity failure.
    # ============================================================

    p <- 75

    mu_t_fun <- function(x) {
      shared_hd_event_location(x)
    }

    sigma_t_fun <- function(x) {

      x <- as.data.frame(x)

      sigma_t <- rep(
        NA_real_,
        nrow(x)
      )

      group_00 <- (
        x$X1 == 0 &
          x$X2 <= 0
      )

      group_01 <- (
        x$X1 == 0 &
          x$X2 > 0
      )

      group_10 <- (
        x$X1 == 1 &
          x$X2 <= 0
      )

      group_11 <- (
        x$X1 == 1 &
          x$X2 > 0
      )

      sigma_t[group_00] <- 0.30
      sigma_t[group_01] <- 0.45
      sigma_t[group_10] <- 0.55
      sigma_t[group_11] <- 0.85

      if (anyNA(sigma_t)) {
        stop(
          paste0(
            "Some observations were not assigned ",
            "an event scale."
          )
        )
      }

      sigma_t
    }

    gen_t <- function(x) {

      draw_min_extreme_value_time(
        mu_t = mu_t_fun(x),
        sigma_t = sigma_t_fun(x)
      )
    }

    gen_c <- function(x) {

      mu_c <- mild_intersection_censoring_location(x)

      exp(
        mu_c +
          0.50 * rnorm(nrow(x))
      )
    }
  } else if (
    setting == "intersection_early_event_mixture_hd"
  ) {

    # ============================================================
    # Shared high-dimensional event location with an early-event
    # mixture whose probability differs across intersections.
    #
    # Early-event probabilities:
    #
    # X1 = 0, X2 <= 0: 0.03
    # X1 = 0, X2 >  0: 0.07
    # X1 = 1, X2 <= 0: 0.10
    # X1 = 1, X2 >  0: 0.18
    #
    # Both components retain the same shared high-dimensional
    # location. The early component is shifted downward.
    # ============================================================

    p <- 75

    mu_t_fun <- function(x) {
      shared_hd_event_location(x)
    }

    early_probability_fun <- function(x) {

      x <- as.data.frame(x)

      prob_early <- rep(
        NA_real_,
        nrow(x)
      )

      group_00 <- (
        x$X1 == 0 &
          x$X2 <= 0
      )

      group_01 <- (
        x$X1 == 0 &
          x$X2 > 0
      )

      group_10 <- (
        x$X1 == 1 &
          x$X2 <= 0
      )

      group_11 <- (
        x$X1 == 1 &
          x$X2 > 0
      )

      prob_early[group_00] <- 0.03
      prob_early[group_01] <- 0.07
      prob_early[group_10] <- 0.10
      prob_early[group_11] <- 0.18

      if (anyNA(prob_early)) {
        stop(
          paste0(
            "Some observations were not assigned ",
            "an early-event probability."
          )
        )
      }

      prob_early
    }

    gen_t <- function(x) {

      n <- nrow(x)

      mu_t <- mu_t_fun(x)

      prob_early <- early_probability_fun(x)

      early_component <- rbinom(
        n = n,
        size = 1,
        prob = prob_early
      )

      u <- pmin(
        pmax(
          runif(n),
          1e-12
        ),
        1 - 1e-12
      )

      eps_t <- log(
        -log(u)
      )

      # The early component has both an earlier location and a
      # slightly narrower distribution.
      component_scale <- ifelse(
        early_component == 1,
        0.25,
        0.35
      )

      early_shift <- 1.35

      log_t <- (
        mu_t -
          early_shift * early_component +
          component_scale * eps_t
      )

      exp(log_t)
    }

    gen_c <- function(x) {

      mu_c <- mild_intersection_censoring_location(x)

      exp(
        mu_c +
          0.50 * rnorm(nrow(x))
      )
    }
  } else if (
    setting == "intersection_location_shift_ld"
  ) {

    # ============================================================
    # Low-dimensional setting with a shared event model and
    # intersection-specific location shifts.
    #
    # Dimension: p = 5
    #
    # The common event slopes are shared across all four groups,
    # but the event location differs according to:
    #
    #   X1 = 0, X2 <= 0:  delta =  0.30
    #   X1 = 0, X2 >  0:  delta = -0.20
    #   X1 = 1, X2 <= 0:  delta =  0.55
    #   X1 = 1, X2 >  0:  delta = -0.65
    #
    # Censoring is mild and does not depend directly on X1.
    # ============================================================

    p <- 5

    intersection_shift_fun <- function(x) {

      x <- as.data.frame(x)

      delta_r <- rep(
        NA_real_,
        nrow(x)
      )

      group_00 <- (
        x$X1 == 0 &
          x$X2 <= 0
      )

      group_01 <- (
        x$X1 == 0 &
          x$X2 > 0
      )

      group_10 <- (
        x$X1 == 1 &
          x$X2 <= 0
      )

      group_11 <- (
        x$X1 == 1 &
          x$X2 > 0
      )

      delta_r[group_00] <- 0.30
      delta_r[group_01] <- -0.20
      delta_r[group_10] <- 0.55
      delta_r[group_11] <- -0.65

      if (anyNA(delta_r)) {
        stop(
          paste0(
            "Some observations were not assigned ",
            "an intersectional location shift."
          )
        )
      }

      delta_r
    }

    mu_t_fun <- function(x) {

      x <- as.data.frame(x)

      2.60 +
        1.00 * x$X2 -
        0.80 * x$X3 +
        0.60 * x$X4 -
        0.40 * x$X5 +
        intersection_shift_fun(x)
    }

    gen_t <- function(x) {

      x <- as.data.frame(x)
      n <- nrow(x)

      # Minimum extreme-value error corresponding to
      # survreg(..., dist = "weibull").
      u <- pmin(
        pmax(
          runif(n),
          1e-12
        ),
        1 - 1e-12
      )

      eps_t <- log(
        -log(u)
      )

      sigma_t <- 0.35

      exp(
        mu_t_fun(x) +
          sigma_t * eps_t
      )
    }

    gen_c <- function(x) {

      x <- as.data.frame(x)

      mu_c <- 3.40 +
        0.15 * x$X3 -
        0.10 * as.numeric(
          x$X2 > 0
        )

      exp(
        mu_c +
          0.50 * rnorm(nrow(x))
      )
    }
  } else if (
    setting == "rare_intersection_scale_shift_ld"
  ) {

    # ============================================================
    # Low-dimensional, rare-intersection scale-shift setting
    #
    # Intended configuration:
    #   p = 5
    #   bernoulli_prob = 0.10
    #
    # The event-location model is shared across all four groups.
    # Only the Weibull scale varies across intersectional groups:
    #
    #   X1 = 0, X2 <= 0: sigma = 0.30
    #   X1 = 0, X2 >  0: sigma = 0.45
    #   X1 = 1, X2 <= 0: sigma = 0.55
    #   X1 = 1, X2 >  0: sigma = 0.80
    #
    # Censoring is moderate and simple. The purpose is to test
    # whether CAMS benefits from pooling the shared event-location
    # structure when minority fitting groups are small.
    # ============================================================

    p <- 5

    mu_t_fun <- function(x) {

      x <- as.data.frame(x)

      required_names <- paste0(
        "X",
        1:5
      )

      missing_names <- setdiff(
        required_names,
        colnames(x)
      )

      if (length(missing_names) > 0L) {
        stop(
          sprintf(
            "Event model is missing columns: %s",
            paste(
              missing_names,
              collapse = ", "
            )
          )
        )
      }

      2.60 +
        0.80 * x$X1 +
        0.80 * x$X2 -
        0.60 * x$X3 +
        0.40 * x$X4 -
        0.30 * x$X5
    }

    sigma_t_fun <- function(x) {

      x <- as.data.frame(x)

      sigma_t <- rep(
        NA_real_,
        nrow(x)
      )

      group_00 <- (
        x$X1 == 0 &
          x$X2 <= 0
      )

      group_01 <- (
        x$X1 == 0 &
          x$X2 > 0
      )

      group_10 <- (
        x$X1 == 1 &
          x$X2 <= 0
      )

      group_11 <- (
        x$X1 == 1 &
          x$X2 > 0
      )

      sigma_t[group_00] <- 0.30
      sigma_t[group_01] <- 0.45
      sigma_t[group_10] <- 0.55
      sigma_t[group_11] <- 0.80

      if (anyNA(sigma_t)) {
        stop(
          paste0(
            "Some observations were not assigned ",
            "an intersection-specific event scale."
          )
        )
      }

      sigma_t
    }

    gen_t <- function(x) {

      draw_min_extreme_value_time(
        mu_t = mu_t_fun(x),
        sigma_t = sigma_t_fun(x)
      )
    }

    gen_c <- function(x) {

      x <- as.data.frame(x)

      mu_c <- (
        3.30 +
          0.10 * x$X3
      )

      exp(
        mu_c +
          0.50 * rnorm(nrow(x))
      )
    }
  } else {
    stop(sprintf("Unknown setting: %s", setting))
  }

  # =====================================================================
  # 2. GENERATE ALL DATA (Unified Logic)
  # =====================================================================
  n_total <- n_train + n_calib + n_test

  # Generate Protected Attribute X1
  X_bern <- rbinom(n_total, 1, bernoulli_prob)

  # Generate Continuous Attributes X2 through Xp
  if (p > 1) {
    X_cont <- matrix(runif(n_total * (p - 1), xmin, xmax), nrow = n_total, ncol = p - 1)
    X <- data.frame(X1 = X_bern, X_cont)
  } else {
    X <- data.frame(X1 = X_bern)
  }

  # Enforce consistent standard column names.
  colnames(X) <- paste0(
    "X",
    1:p
  )

  # Construct the supplied nonlinear basis for this setting.
  if (
    setting == "basis_intersection_shared_weibull"
  ) {

    X$X4 <- X$X2^2
    X$X5 <- sin(X$X3)
  }

  # Calculate event and censoring times.
  T_time <- gen_t(X)
  C_time <- gen_c(X)

  event <- (T_time < C_time)
  censored_T <- pmin(T_time, C_time)

  # Bind everything into a master data frame
  data_full <- data.frame(X, C = C_time, censored_T = censored_T, event = event)

  # =====================================================================
  # 3. SPLIT DATA INTO TRAIN, CALIB, AND TEST
  # =====================================================================
  data_fit <- data_full[1:n_train, , drop=FALSE]
  data_calib <- data_full[(n_train + 1):(n_train + n_calib), , drop=FALSE]
  data_test <- data_full[(n_train + n_calib + 1):n_total, , drop=FALSE]

  # Original code expects 'data' to contain fit + calib 

  # Extract exact test event times
  T_test <- T_time[(n_train + n_calib + 1):n_total]

  # Collect results 
  obj <- list(data_fit = data_fit,
              data_calib = data_calib,
              data_test = data_test,
              T_test = T_test,
              p = p)

  return(obj)
}
