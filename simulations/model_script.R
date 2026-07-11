model_generating_fun <- function(n_train, n_calib, n_test,
                                 setting, xmin, xmax,
                                 bernoulli_prob = 0.1) {

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

  # Enforce consistent standard column names (X1, X2, ... Xp)
  colnames(X) <- paste0("X", 1:p)

  # Calculate Event Times and Censoring Times based on the selected setting formulas
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