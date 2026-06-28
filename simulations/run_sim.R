total_start_time <- proc.time()[3]
########################################
## Process the input argument
########################################
args <- commandArgs(trailingOnly = TRUE)
setting_list <- unlist(strsplit(args[1], ","))
# setting_list <- c("starve_hetero_high_dim")
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

alpha <- 0.1    # target level 1-alpha
n <- 1000
n_test <- 5000
n_train <- n
n_calib <- n
xmin <- -2
xmax <- 2
bernoulli_prob <- 0.1

num_runs <- 5

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
    run_folder <- sprintf("../new_results%s/%d", as.character(bernoulli_prob), current_seed)
  } else {
    run_folder <- sprintf("../results%s/%d", as.character(bernoulli_prob), current_seed)
  }
  dir.create(run_folder, showWarnings = FALSE, recursive = TRUE)
  run_start_time <- proc.time()[3]
  j <- 1
  for(setting in setting_list){
    cat(sprintf("\n=== Run %d | Setting: %s ===\n", i, setting))

    start_time <- proc.time()[3]
    # Run the simulation
    simures <- simu(current_seed + 1234, setting, only_cams,
                    n_train, n_calib, n_test,
                    xmin, xmax, alpha,
                    bernoulli_prob)
    cat(sprintf("%s for run %d: %.2f seconds.\n", setting, i, proc.time()[3] - start_time))
    # Save the result file directly into the newly created numbered folder
    save_dir <- sprintf("%s/%d. %s_seed_%d.csv", run_folder, j, setting, current_seed)
    write.csv(simures, save_dir)
    j <- j + 1
  }

  cat(sprintf("\nCompleted run %d/%d in %.2f seconds.\n", i, num_runs, proc.time()[3] - run_start_time))
}

if (only_cams) {
  merge()
}

plots_dir <- sprintf("../plots%s", as.character(bernoulli_prob))
make_plots(results_dir = run_folder,
           plots_dir = plots_dir,
           target_alpha = alpha)

cat(sprintf("\nCompleted in %.2f seconds.\n", proc.time()[3] - total_start_time))