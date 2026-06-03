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
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(snow)) # Needed for detecting cores

########################################
### source code
########################################
source("./source_code.R")
source("./model_script.R")
source("./simu.R")

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
n <- 500
n_test <- 2500
n_train <- n
n_calib <- n
xmin <- 0 
xmax <- 4
beta <- 20 / sqrt(n)
exp_rate <- .1

num_runs <- 5

# Detect cores just to print a helpful message (the actual multithreading happens inside the utils scripts)
slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
num_cores <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
cat(sprintf("Starting sequential outer loop. Inner algorithms will utilize %d cores...\n", num_cores))

########################################
### SEQUENTIAL LOOP (Letting inner functions multithread)
########################################
for(i in 1:num_runs){
  
  # Update the seed for each run
  current_seed <- i 
  
  # Create a separate folder named with the run number inside the 'results' folder
  run_folder <- sprintf("../results/%d", current_seed)
  dir.create(run_folder, showWarnings = FALSE, recursive = TRUE)
  
  for(setting in setting_list){
    cat(sprintf("\n=== Run %d | Setting: %s ===\n", i, setting))
    
    if(setting %in% c("hd_homosc","hd_heterosc")){
      p <- 10
    }else{
      p <- 1
    }
    
    # Run the simulation
    simures <- simu(current_seed + 1234, setting, n, p,
                    n_train, n_calib, n_test,
                    beta, xmin, xmax,
                    exp_rate, alpha)
    
    # Save the result file directly into the newly created numbered folder
    save_dir <- sprintf("%s/%s_seed_%d.csv", run_folder, setting, current_seed)
    write.csv(simures, save_dir)
  }
  
  cat(sprintf("\nCompleted run %d/%d\n", i, num_runs))
}