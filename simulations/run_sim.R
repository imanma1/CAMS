########################################
## Process the input argument
########################################
args <- commandArgs(trailingOnly = TRUE)
seed<- as.integer(args[1])
if(is.na(seed)){seed <- 1}

########################################
## load libraries
########################################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(quantreg))
suppressPackageStartupMessages(library(GauPro))
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(grf))

# NEW: Load Parallel Libraries
library(doParallel)
library(foreach)

########################################
### run simulations
########################################
## configurations
setting_list = c("ld_setting1",
                 "ld_setting2",
                 "ld_setting3",
                 "ld_setting4",
                 "hd_homosc",
                 "hd_heterosc")

alpha <- .1    # target level 1-alpha
n <- 1200
n_test <- 6000
n_train <- n
n_calib <- n
xmin <- 0 
xmax <- 4
beta <- 20 / sqrt(n)
exp_rate <- .1

# Define the number of experimental runs
num_runs <- 5

########################################
### SETUP PARALLEL CLUSTER
########################################
# Use all available cores minus 1 to keep your computer responsive
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat(sprintf("Starting parallel simulation across %d cores...\n", num_cores))

########################################
### PARALLEL LOOP
########################################
# Replace the standard 'for' loop with 'foreach %dopar%'
foreach(i = 1:num_runs,
        .packages = c("tidyverse", "survival", "quantreg", "GauPro", "gbm", "grf")) %dopar% {
  
  # IMPORTANT: Source the custom scripts inside the loop so each parallel worker loads them
  source("./source_code.R")
  source("./model_script.R")
  source("./simu.R")
  
  # Update the seed for each run
  current_seed <- i 
  
  # Create a separate folder named with the run number inside the 'results' folder
  run_folder <- sprintf("../results/%d", current_seed)
  dir.create(run_folder, showWarnings = FALSE, recursive = TRUE)
  
  for(setting in setting_list){
    if(setting %in% c("hd_homosc","hd_heterosc")){
      p <- 10
    }else{
      p <- 1
    }
    
    # Run the simulation for the current setting and seed
    simures <- simu(current_seed + 1234, setting, n, p,
                    n_train, n_calib, n_test,
                    beta, xmin, xmax,
                    exp_rate, alpha)
    
    # Save the result file directly into the newly created numbered folder
    save_dir <- sprintf("%s/%s_seed_%d.csv", run_folder, setting, current_seed)
    write.csv(simures, save_dir)
  }
}

########################################
### CLEANUP
########################################
stopCluster(cl)
cat("All parallel runs completed successfully!\n")