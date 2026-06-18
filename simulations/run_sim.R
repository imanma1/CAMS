total_start_time <- proc.time()[3]
########################################
## Process the input argument
########################################
args <- commandArgs(trailingOnly = TRUE)
setting_list <- unlist(strsplit(args[1], ","))
seed <- as.integer(args[2])
if (is.na(seed)) {
  seed <- 1
}
only_cams <- as.logical(as.integer(args[3]))
if (is.na(only_cams)) {
  only_cams <- FALSE
}

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
source("./merge.R")
source("./fig.R")

########################################
### run simulations
########################################
## configurations
########################################

# setting_list <- c("starve_hetero_high_dim")

alpha <- .1    # target level 1-alpha
n <- 1000
n_test <- 5000
n_train <- n
n_calib <- n
xmin <- -2
xmax <- 2

num_runs <- 1

# Detect cores just to print a helpful message (the actual multithreading happens inside the utils scripts)
slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
num_cores <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
cat(sprintf("Starting sequential outer loop. Inner algorithms will utilize %d cores...\n", num_cores))

########################################
### SEQUENTIAL LOOP (Letting inner functions multithread)
########################################
for (i in 1:num_runs) {
  # Update the seed for each run
  current_seed <- seed + (i - 1)
  
  # Create a separate folder named with the run number inside the 'results' folder
  if (only_cams){
    run_folder <- sprintf("../new_results/%d", current_seed)
  } else {
     run_folder <- sprintf("../results/%d", current_seed)
  }
  dir.create(run_folder, showWarnings = FALSE, recursive = TRUE)
  run_start_time <- proc.time()[3]
  for(setting in setting_list){
    cat(sprintf("\n=== Run %d | Setting: %s ===\n", i, setting))
    
    start_time <- proc.time()[3]
    # Run the simulation
    simures <- simu(current_seed + 1234, setting, only_cams,
                    n_train, n_calib, n_test,
                    xmin, xmax, alpha)
    cat(sprintf("%s for run %d: %.2f seconds.\n", setting, i, proc.time()[3] - start_time))
    
    # Save the result file directly into the newly created numbered folder
    save_dir <- sprintf("%s/%s_seed_%d.csv", run_folder, setting, current_seed)
    write.csv(simures, save_dir)
  }
  
  cat(sprintf("\nCompleted run %d/%d in %.2f seconds.\n", i, num_runs, proc.time()[3] - run_start_time))
}

if (only_cams) {
  merge()
}
make_plots(setting_list)

cat(sprintf("\nCompleted in %.2f seconds.\n", proc.time()[3] - total_start_time))