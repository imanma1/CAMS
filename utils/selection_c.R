#' Automatic selection of censoring time
#'
#' The funciton that automatically choose a value of c_0.
#'
#' @export

selection_c <- function(X,C,event,time,alpha,
                        c_ref,weight_ref,
                        model="cox",
                        type="quantile",
                        dist="weibull",
                        mdl0=NULL){ 
  
  ## Get the dimension of the input
  if(is.null(dim(X))){
    n <- length(X)
    p <- 1
  }else{
    n <- dim(X)[1]
    p <- dim(X)[2]
  }
  xnames <- paste0("X",1:p)
  data <- cbind(X,C,event,time)
  data <- data.frame(data)
  colnames(data) <- c(xnames,"C","event","censored_T")

  # === OPTIMIZATION: PARALLELIZE GRID SEARCH ===
  library(parallel)
  if (.Platform$OS.type == "windows") {
    n_threads <- 1  
  } else {
    slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
    n_threads <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
  }

  ## Evaluate the average bound for each candidate c simultaneously
  bnd_ref <- unlist(mclapply(1:length(c_ref), function(i) {
    
    # Safely handle the weight_ref subsetting
    current_weight <- if(is.null(weight_ref)) NULL else weight_ref[,i]
    
    # Calculate and return the bound for this specific threshold
    evaluate_length(c_ref[i], alpha=alpha, n=n, p=p, model=model,
                    data=data, weight=current_weight, xnames=xnames,
                    type=type, dist=dist, mdl0=mdl0)
                    
  }, mc.cores = n_threads))

  # Find the optimal cut-off based on the max bound
  c_opt <- c_ref[which.max(bnd_ref)]
  return(list(c_opt=c_opt, c_ref=c_ref, bnd_ref=bnd_ref))
}

evaluate_length <- function(c,alpha,n,p,
                            model, data, weight, xnames,
                            type = "quantile",
                            dist = "weibull",
                            seed = 2020,
                            mdl0 = NULL){ # <--- ADD mdl0
  
  set.seed(seed)
  I_fit <- sample(1:n,floor(n/2),replace=FALSE)
  I_calib <- sample((1:n)[-I_fit],floor(n/4),replace=FALSE)
  I_test <- (1:n)[-c(I_fit,I_calib)]
  
  data_fit <- data[I_fit,]
  data_calib <- data[I_calib,]
  data_test <- data[I_test,]
  
  if(is.null(weight)){
    # === OPTIMIZATION 3: REUSE mdl0 ===
    # Use cens_prob (which uses pre-trained mdl0) instead of censoring_prob (which trains from scratch)
    res <- cens_prob(mdl=mdl0, calib=data_calib, test=data_test,
                     method="gpr", xnames=xnames, c=c)
    pr_calib <- res$pr_calib
    pr_new <- res$pr_new
    weight_calib <- 1/pr_calib
    weight_new <- 1/pr_new
  }else{
    weight_calib <- weight[I_calib]
    weight_new <- weight[I_test]
  }
  x <- data_test[,colnames(data_test)%in%xnames, drop=FALSE]
  
  if(model == "cox"){
    bnd <- cox0_based(x,c,alpha,
                     data_fit,
                     data_calib,
                     type = "quantile",
                     dist,
                     weight_calib,
                     weight_new)
   }
  
  if(model == "randomforest"){
    bnd <- rf0_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new)
  }
  
  if(model == "pow"){
    bnd <- pow0_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new)
  }

  if(model == "portnoy"){
    bnd <- portnoy0_based(x,c,alpha,
                        data_fit,
                        data_calib,
                        weight_calib,
                        weight_new)
  }

  if(model == "PengHuang"){
    bnd <- ph0_based(x,c,alpha,
                   data_fit,
                   data_calib,
                   weight_calib,
                   weight_new)
  }

  return(mean(bnd))
}
