# Added a bernouli(0.3) variable for all settings as the first variable

model_generating_fun <- function(n_train, n_calib, n_test,
                                 setting, beta, xnames, xmin, xmax, exp_rate){

  bernoulli_prob <- 0.3
  
  if(setting == "ld_setting1"){
    p <- 2 # Updated to 2 features
    sigma_x <- function(x) (5 - x[,2])/10 # Shifted to x[,2]
    
    ########################################
    ## Data generating models
    ########################################
    gen_t <- function(x) exp(beta * x[,2] +  2 * rnorm(nrow(x)))
    gen_c <- function(x) rexp(rate = exp_rate, n = nrow(x))
    
    ########################################
    ## Generate training data
    ########################################
    X_bern <- rbinom(n_train, 1, bernoulli_prob)
    X_cont <- runif(n_train, xmin, xmax)
    X <- data.frame(X1 = X_bern, X2 = X_cont)
    
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data_fit <- data.frame(X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X_bern <- rbinom(n_calib + n_test, 1, bernoulli_prob)
    X_cont <- runif(n_calib + n_test, xmin, xmax)     
    X <- data.frame(X1 = X_bern, X2 = X_cont)
    
    T <- gen_t(X) 
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X, C = C, event = event, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  
  if(setting == "ld_setting2"){
    p <- 2
    
    ########################################
    ## Data generating models
    ########################################
    sigma_x <- function(x) (5 + x[,2])/10
    gen_t <- function(x) exp(3*(x[,2]>2) + 1*x[,2]*(x[,2]<=2) + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) rexp(rate = exp_rate, n = nrow(x))
    
    ########################################
    ## Generate training data
    ########################################
    X_bern <- rbinom(n_train, 1, bernoulli_prob)
    X_cont <- runif(n_train, xmin, xmax)
    X <- data.frame(X1 = X_bern, X2 = X_cont)
    
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data_fit <- data.frame(X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X_bern <- rbinom(n_calib + n_test, 1, bernoulli_prob)
    X_cont <- runif(n_calib + n_test, xmin, xmax)     
    X <- data.frame(X1 = X_bern, X2 = X_cont)
    
    T <- gen_t(X) 
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X, C = C, event = event, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  
  if(setting == "ld_setting3"){
    p <- 2
    ########################################
    ## Data generating models
    ########################################
    gen_t <- function(x) exp(2 * (x[,2]>2) + 1 * x[,2] *(x[,2]<=2) + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) rexp(rate = exp_rate * (2.5 + (6+x[,2])/10), n = nrow(x))
    
    ########################################
    ## Generate training data
    ########################################
    X_bern <- rbinom(n_train, 1, bernoulli_prob)
    X_cont <- runif(n_train, xmin, xmax)
    X <- data.frame(X1 = X_bern, X2 = X_cont)
    
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data_fit <- data.frame(X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X_bern <- rbinom(n_calib + n_test, 1, bernoulli_prob)
    X_cont <- runif(n_calib + n_test, xmin, xmax)     
    X <- data.frame(X1 = X_bern, X2 = X_cont)
    
    T <- gen_t(X) 
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X, C = C, event = event, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  
  if(setting == "ld_setting4"){
    p <- 2
    
    ########################################
    ## Data generating models
    ########################################
    gen_t <- function(x) exp(3 * (x[,2]>2) + 1.5 * x[,2] *(x[,2]<=2) + 0.5 * rnorm(nrow(x)))
    gen_c <- function(x) exp(2 + (2-x[,2]) / 50 + 0.5 * rnorm(nrow(x)))
    
    ########################################
    ## Generate training data
    ########################################
    X_bern <- rbinom(n_train, 1, bernoulli_prob)
    X_cont <- runif(n_train, xmin, xmax)
    X <- data.frame(X1 = X_bern, X2 = X_cont)
    
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data_fit <- data.frame(X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X_bern <- rbinom(n_calib + n_test, 1, bernoulli_prob)
    X_cont <- runif(n_calib + n_test, xmin, xmax)     
    X <- data.frame(X1 = X_bern, X2 = X_cont)
    
    T <- gen_t(X) 
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X, C = C, event = event, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  
  if(setting == "hd_homosc"){
    p_cont <- length(xnames)
    new_xnames <- paste0("X", 1:(p_cont + 1)) # Create updated names for p+1 variables
    
    ########################################
    ## Data generating models (Indices shifted +1)
    ########################################
    mu_x <- function(x) (beta * x[,2] + beta * sqrt(x[,4] * x[,6])) / 5 + 1
    gen_t <- function(x) exp(mu_x(x) + rnorm(nrow(x)))
    gen_c <- function(x) rexp(rate = exp_rate * (x[,11] + 0.5), n = nrow(x))
    
    ## Generate training data
    X_cont <- matrix(runif(n_train * p_cont, min = xmin, max = xmax), n_train)
    X_bern <- rbinom(n_train, 1, bernoulli_prob)
    X <- cbind(X_bern, X_cont)
    
    T <- gen_t(X)
    C <- gen_c(X) 
    event <- (T<C)
    censored_T <- pmin(T,C)
    data_fit <- data.frame(X, C = C, censored_T = censored_T, event = event)
    colnames(data_fit) <- c(new_xnames, "C", "censored_T", "event")
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X_cont <- matrix(runif((n_calib + n_test) * p_cont, min = xmin, max = xmax), n_calib + n_test)
    X_bern <- rbinom(n_calib + n_test, 1, bernoulli_prob)
    X <- cbind(X_bern, X_cont)
    
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T<C)
    censored_T <- pmin(T,C)
    data <- data.frame(X, C = C, censored_T = censored_T,  event = event)
    colnames(data) <- c(new_xnames, "C", "censored_T", "event")
    data_calib <- data[1:n_calib,]
    data_test <- data[(n_calib+1) : (n_calib+n_test),]
    data <- rbind(data_fit,data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  
  if(setting == "hd_heterosc"){
    p_cont <- length(xnames)
    new_xnames <- paste0("X", 1:(p_cont + 1))
    
    ########################################
    ## Data generating models (Indices shifted +1)
    ########################################
    mu_x <- function(x) (beta * x[,2] + beta * sqrt(x[,4] * x[,6])) / 5 + 1
    sigma_x <- function(x) (x[,3] + 2) / 4 # Shifted from x[,2] to x[,3]
    gen_t <- function(x) exp(mu_x(x) + sigma_x(x) * rnorm(nrow(x)))
    gen_c <- function(x) rexp(rate = exp_rate * (x[,11] + 0.5), n = nrow(x))
    
    ## Generate training data
    X_cont <- matrix(runif(n_train * p_cont, min = xmin, max = xmax), n_train)
    X_bern <- rbinom(n_train, 1, bernoulli_prob)
    X <- cbind(X_bern, X_cont)
    
    T <- gen_t(X)
    C <- gen_c(X) 
    event <- (T<C)
    censored_T <- pmin(T,C)
    data_fit <- data.frame(X, C = C, censored_T = censored_T, event = event)
    colnames(data_fit) <- c(new_xnames, "C", "censored_T", "event")
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    X_cont <- matrix(runif((n_calib + n_test) * p_cont, min = xmin, max = xmax), n_calib + n_test)
    X_bern <- rbinom(n_calib + n_test, 1, bernoulli_prob)
    X <- cbind(X_bern, X_cont)
    
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T<C)
    censored_T <- pmin(T,C)
    data <- data.frame(X, C = C, censored_T = censored_T,  event = event)
    colnames(data) <- c(new_xnames, "C", "censored_T", "event")
    data_calib <- data[1:n_calib,]
    data_test <- data[(n_calib+1) : (n_calib+n_test),]
    data <- rbind(data_fit,data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  
  ## Collect results 
  obj <- list(data_fit = data_fit, 
              data_calib = data_calib, 
              data_test = data_test, 
              data = data, 
              T_test = T_test)
  
  return(obj)
}