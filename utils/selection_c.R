selection_c <- function(data, p, n, xnames,
                        c_ref, weight_ref,
                        model = "cox",
                        type = "quantile",
                        dist = "weibull",
                        mdl0 = NULL,
                        alpha,
                        use_oracle_sc = FALSE) {
  # === OPTIMIZATION: PARALLELIZE GRID SEARCH ===
  library(parallel)
  if (.Platform$OS.type == "windows") {
    n_threads <- 1
  } else {
    slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    n_threads <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
  }
  
  if (.Platform$OS.type != "windows") {
    library(RhpcBLASctl)
    original_threads <- blas_get_num_procs() 
    blas_set_num_threads(1)
    omp_set_num_threads(1)
  } 

  ## Evaluate the average bound for each candidate c simultaneously
  bnd_ref <- unlist(mclapply(1:length(c_ref), function(i) {
    
    current_weight <- if(is.null(weight_ref)) NULL else weight_ref[,i]
    
    evaluate_length(c_ref[i],
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
                    use_oracle_sc = use_oracle_sc)
                    
  }, mc.cores = n_threads))

  if (.Platform$OS.type != "windows") {
    blas_set_num_threads(original_threads)
    omp_set_num_threads(original_threads)
  }

  c_opt <- c_ref[which.max(bnd_ref)]
  return(list(c_opt = c_opt, c_ref = c_ref, bnd_ref = bnd_ref))
}


evaluate_length <- function(c, alpha, n, p,
                            model, data, weight, xnames,
                            type = "quantile",
                            dist = "weibull",
                            seed = 2020,
                            mdl0 = NULL,
                            use_oracle_sc = FALSE) {
  
  set.seed(seed)

  I_fit <- sample(1:n, floor(n/2), replace = FALSE)
  I_calib <- sample((1:n)[-I_fit], floor(n/4), replace = FALSE)
  I_test <- (1:n)[-c(I_fit, I_calib)]
  
  data_fit <- data[I_fit, ]
  data_calib <- data[I_calib, ]
  data_test <- data[I_test, ]
  
  if (is.null(weight)) {

    if (use_oracle_sc || inherits(mdl0, "oracle_sc")) {

      pr_calib <- sc_prob(
        mdl0 = mdl0,
        data = data_calib,
        xnames = xnames,
        t = c
      )

      pr_new <- sc_prob(
        mdl0 = mdl0,
        data = data_test,
        xnames = xnames,
        t = c
      )

    } else {

      res <- cens_prob(mdl = mdl0,
                       calib = data_calib,
                       test = data_test,
                       method = "gpr",
                       xnames = xnames,
                       c = c)

      pr_calib <- res$pr_calib
      pr_new <- res$pr_new
    }

    weight_calib <- 1 / pr_calib
    weight_new <- 1 / pr_new

  } else {
    weight_calib <- weight[I_calib]
    weight_new <- weight[I_test]
  }

  x <- data_test[, colnames(data_test) %in% xnames, drop = FALSE]

  bnd <- cox0_based(
    x = x,
    p = p,
    len_x = nrow(x),
    xnames = xnames,
    c = c,
    alpha = alpha,
    data_fit = data_fit,
    data_calib = data_calib,
    type = "quantile",
    dist = dist,
    weight_calib = weight_calib,
    weight_new = weight_new
  )

  return(mean(bnd))
}