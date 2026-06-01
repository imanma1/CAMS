# cams <- function(x,
#                      Xtrain, C, event, time,
#                      alpha=0.1, p,
#                      mdl0
# ){

#   if(is.null(dim(x)[1])){
#     len_x <- length(x)
#     p <- 1
#   }else{
#     len_x <- dim(x)[1]
#     p <- dim(x)[2]
#   }

#   X <- Xtrain
#   xnames <- paste0('X', 1:p)
#   data <- as.data.frame(cbind(C,event,time,X))
#   colnames(data) <- c("C","event","censored_T",xnames)

#   ## Split the data into the training set and the calibration set
#   n <- dim(data)[1]
#   n_train <- n/2
#   n_calib <- n - n_train
#   I_fit <- sample(1:n, n_train, replace = FALSE)
#   data_fit <- data[I_fit,]
#   data_calib <- data[-I_fit,]

#   n <- nrow(data_calib)
#   ## Fit the survival model
#   xnames <- paste0("X",1:p)
#   newdata <- data.frame(x)
#   colnames(newdata) <- xnames
#   fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
#   mdl <- survreg(fmla, data = data_fit, dist= "weibull")
  
#   # cutoff for quantile of C
#   # to bound weights
#   cens_rt = 1-1/log(n)
  
#   start_time = proc.time()[3]
#   qt_res <- alpha_qt(mdl, newdata,
#                      data_fit, data_calib,
#                      xnames, alpha, len_x,
#                      cens_rt = cens_rt,
#                      mdl0 = mdl0)
#   lower_bnd_qtg <- qt_res$lower_bnd_g
#   lower_bnd_qtl <- qt_res$lower_bnd_l
#   end_time <- proc.time()[3]
#   time_qt <- end_time - start_time
# }


# alpha_qt <- function(mdl, newdata,
#                      data_fit, data_calib,
#                      xnames, alpha, len_x,
#                      mdl0, cens_rt){
  
#   v_list = v_pts_qt(mdl,
#                       data_fit, data_calib,
#                       xnames, alpha, cens_rt)
#   v_list = sort(unique(as.numeric(v_list)))
  
#   ## Obtain the final confidence interval
#   lower_bnd_l <- rep(0,len_x)
#   lower_bnd_g <- rep(0,len_x)
  
#   alpha_v_list <- sapply(v_list, est_alpha_qt, mdl = mdl, 
#                     data_calib = data_calib, xnames = xnames, alpha = alpha, 
#                     cens_rt = cens_rt, mdl0 = mdl0, newdata = newdata)

#   # monotonize alpha
#   alpha_v <- monot(alpha_v_list)
#   # return 0 if alpha is above target level
#   if(sum(alpha_v<=alpha)==0){
#     v_hat_l = NULL
#     v_hat_g = NULL
#   }else{
#     v_hat_l <- min(v_list[alpha_v <= alpha])
#   }
  
#   for(i in 1:len_x){
#     nxi <- data.frame(newdata[i,])
#     colnames(nxi) <- xnames
    
#     lower_bnd_l[i] <- lv_qt(mdl, nxi, v_hat_l, alpha, cens_rt)  
#   }
  
#   return(list(lower_bnd_l = lower_bnd_l, 
#               lower_bnd_g = lower_bnd_g))
  
# }