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
########################################
setting_list = c("homo_cens", "cov_cens", "prot_cens", 
                 "heavy_prot_cens", "heavy_inter_cens", 
                 "surv_misspec", "cens_misspec", "simul_misspec")

alpha <- .1    # target level 1-alpha
n <- 1000
n_test <- 5000
n_train <- n
n_calib <- n
xmin <- -2
xmax <- 2

num_runs <- 50

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
    
    # Run the simulation
    simures <- simu(current_seed + 1234, setting,
                    n_train, n_calib, n_test,
                    xmin, xmax, alpha)
    
    # Save the result file directly into the newly created numbered folder
    save_dir <- sprintf("%s/%s_seed_%d.csv", run_folder, setting, current_seed)
    write.csv(simures, save_dir)
  }
  
  cat(sprintf("\nCompleted run %d/%d\n", i, num_runs))
}