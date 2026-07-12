total_start_time <- proc.time()[3]
########################################
## Process the input argument
########################################
args <- commandArgs(trailingOnly = TRUE)
default_settings <- c(
  "cams_vs_vanilla_lower_tail_hd_mild",
  "cams_vs_vanilla_lower_tail_hd_main",
  "cams_vs_vanilla_lower_tail_hd_strong"
)
setting_list <- if (length(args) >= 1L && nzchar(args[1])) {
  unlist(strsplit(args[1], ",", fixed = TRUE))
} else {
  default_settings
}

seed <- if (length(args) >= 2L) as.integer(args[2]) else NA_integer_
if (is.na(seed)) {
  seed <- 1
}

sc_arg <- if (length(args) >= 3L) args[3] else "oracle"
if (sc_arg %in% c("0", "1")) {
  # Backward compatibility with the previous use_oracle_sc argument.
  sc_methods <- if (sc_arg == "1") "oracle" else "aft_lognormal"
} else {
  sc_methods <- unlist(strsplit(sc_arg, ",", fixed = TRUE))
}
valid_sc_methods <- c("oracle", "rsf", "aft_lognormal", "km_x1", "km")
if (any(!sc_methods %in% valid_sc_methods)) {
  stop(sprintf(
    "Unknown censoring model(s): %s",
    paste(setdiff(sc_methods, valid_sc_methods), collapse = ", ")
  ))
}

only_cams <- if (length(args) >= 4L) as.logical(as.integer(args[4])) else NA
if (is.na(only_cams)) {
  only_cams <- FALSE
}

augmentation_arg <- if (length(args) >= 5L) args[5] else "correct"
augmentation_methods <- unlist(strsplit(augmentation_arg, ",", fixed = TRUE))
valid_augmentation_methods <- c("same", "correct", "wrong")
if (any(!augmentation_methods %in% valid_augmentation_methods)) {
  stop(sprintf(
    "Unknown augmentation model(s): %s",
    paste(setdiff(augmentation_methods, valid_augmentation_methods), collapse = ", ")
  ))
}

homoscedastic_event <- if (length(args) >= 6L) {
  as.logical(as.integer(args[6]))
} else {
  FALSE
}
if (is.na(homoscedastic_event)) homoscedastic_event <- FALSE

sc_ntree <- if (length(args) >= 7L) as.integer(args[7]) else 1000L
if (is.na(sc_ntree) || sc_ntree < 1L) sc_ntree <- 1000L

num_runs <- if (length(args) >= 8L) as.integer(args[8]) else 10L
if (is.na(num_runs) || num_runs < 1L) num_runs <- 10L

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

# Detect cores just to print a helpful message (the actual multithreading happens inside the utils scripts)
slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
num_cores <- ifelse(is.na(slurm_cores), detectCores(), slurm_cores)
cat(sprintf("Starting sequential outer loop. Inner algorithms will utilize %d cores...\n", num_cores))

########################################
### SEQUENTIAL LOOP (Letting inner functions multithread)
########################################
base_results_dir <- if (only_cams) {
  sprintf("../new_results%s", as.character(bernoulli_prob))
} else {
  sprintf("../results%s", as.character(bernoulli_prob))
}

event_scale_label <- if (homoscedastic_event) "event_constant_scale" else "event_group_varying_scale"

for (i in seq_len(num_runs)) {
  # Update the seed for each run
  current_seed <- seed + (i - 1)

  run_start_time <- proc.time()[3]
  for (sc_method in sc_methods) {
    for (augmentation_method in augmentation_methods) {
      run_folder <- file.path(
        base_results_dir,
        paste0("sc_", sc_method),
        paste0("aug_", augmentation_method),
        event_scale_label,
        as.character(current_seed)
      )
      dir.create(run_folder, showWarnings = FALSE, recursive = TRUE)

      for (j in seq_along(setting_list)) {
        setting <- setting_list[j]
        cat(sprintf(
          "\n=== Run %d | Setting: %s | G: %s | augmentation: %s ===\n",
          i, setting, sc_method, augmentation_method
        ))

        start_time <- proc.time()[3]
        simures <- simu(
          current_seed + 1234,
          setting,
          only_cams,
          n_train,
          n_calib,
          n_test,
          xmin,
          xmax,
          alpha,
          bernoulli_prob,
          use_oracle_sc = identical(sc_method, "oracle"),
          sc_method = sc_method,
          augmentation_method = augmentation_method,
          homoscedastic_event = homoscedastic_event,
          sc_ntree = sc_ntree
        )
        cat(sprintf(
          "%s for run %d (%s/%s): %.2f seconds.\n",
          setting, i, sc_method, augmentation_method,
          proc.time()[3] - start_time
        ))
        save_dir <- file.path(
          run_folder,
          sprintf("%d. %s_seed_%d.csv", j, setting, current_seed)
        )
        write.csv(simures, save_dir, row.names = FALSE)
      }
    }
  }

  cat(sprintf("\nCompleted run %d/%d in %.2f seconds.\n", i, num_runs, proc.time()[3] - run_start_time))
}

for (sc_method in sc_methods) {
  for (augmentation_method in augmentation_methods) {
    combination_results_dir <- file.path(
      base_results_dir,
      paste0("sc_", sc_method),
      paste0("aug_", augmentation_method),
      event_scale_label
    )
    combination_plots_dir <- file.path(
      sub("results", "plots", base_results_dir, fixed = TRUE),
      paste0("sc_", sc_method),
      paste0("aug_", augmentation_method),
      event_scale_label
    )
    make_plots(
      results_dir = combination_results_dir,
      plots_dir = combination_plots_dir,
      target_alpha = alpha
    )
  }
}

cat(sprintf("\nCompleted in %.2f seconds.\n", proc.time()[3] - total_start_time))
