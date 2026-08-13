options(debug_dft_fixed = FALSE)

total_start_time <- proc.time()[3]

########################################
## Process input arguments
########################################

args <- commandArgs(trailingOnly = TRUE)

# First argument format:
#
# setting:n_train:n_calib:n_test:bernoulli_prob
#
# Multiple settings:
#
# setting1:n_train:n_calib:n_test:prob;
# setting2:n_train:n_calib:n_test:prob;
# ...

default_setting_configs <- paste(
  "cams_vs_vanilla_lower_tail_hd_main",
  2000,
  8000,
  20000,
  0.10,
  sep = ":"
)


########################################
## Parse per-setting configurations
########################################

setting_config_arg <- if (
  length(args) >= 1L &&
  nzchar(args[1])
) {
  args[1]
} else {
  default_setting_configs
}


parse_setting_configs <- function(config_string) {

  entries <- unlist(
    strsplit(
      config_string,
      ";",
      fixed = TRUE
    )
  )

  entries <- trimws(entries)
  entries <- entries[nzchar(entries)]

  if (length(entries) == 0L) {
    stop("No setting configurations were provided.")
  }


  configs <- lapply(
    entries,
    function(entry) {

      parts <- trimws(
        unlist(
          strsplit(
            entry,
            ":",
            fixed = TRUE
          )
        )
      )

      if (length(parts) != 5L) {

        stop(
          paste0(
            "Invalid setting configuration: ",
            entry,
            "\nExpected format: ",
            "setting:n_train:n_calib:n_test:bernoulli_prob"
          )
        )
      }


      data.frame(
        setting = parts[1],
        n_train = as.integer(parts[2]),
        n_calib = as.integer(parts[3]),
        n_test = as.integer(parts[4]),
        bernoulli_prob = as.numeric(parts[5]),
        stringsAsFactors = FALSE
      )
    }
  )


  config_df <- do.call(
    rbind,
    configs
  )

  rownames(config_df) <- NULL


  if (any(!nzchar(config_df$setting))) {
    stop(
      "Each configuration must contain a non-empty setting name."
    )
  }


  if (
    anyNA(config_df$n_train) ||
    anyNA(config_df$n_calib) ||
    anyNA(config_df$n_test)
  ) {
    stop(
      "n_train, n_calib, and n_test must be integers."
    )
  }


  if (
    any(config_df$n_train < 1L) ||
    any(config_df$n_calib < 1L) ||
    any(config_df$n_test < 1L)
  ) {
    stop(
      "n_train, n_calib, and n_test must all be positive."
    )
  }


  if (
    anyNA(config_df$bernoulli_prob) ||
    any(!is.finite(config_df$bernoulli_prob)) ||
    any(config_df$bernoulli_prob <= 0) ||
    any(config_df$bernoulli_prob >= 1)
  ) {
    stop(
      "bernoulli_prob must be strictly between 0 and 1."
    )
  }


  if (anyDuplicated(config_df$setting)) {

    duplicated_settings <- unique(
      config_df$setting[
        duplicated(config_df$setting)
      ]
    )

    stop(
      sprintf(
        "Each setting should appear only once. Duplicates: %s",
        paste(
          duplicated_settings,
          collapse = ", "
        )
      )
    )
  }


  config_df
}


setting_configs <- parse_setting_configs(
  setting_config_arg
)


########################################
## Other command-line arguments
########################################

seed <- if (length(args) >= 2L) {
  as.integer(args[2])
} else {
  NA_integer_
}

if (is.na(seed)) {
  seed <- 1L
}


sc_arg <- if (length(args) >= 3L) {
  args[3]
} else {
  "oracle"
}


if (sc_arg %in% c("0", "1")) {

  # Backward compatibility with previous
  # use_oracle_sc argument.

  sc_methods <- if (sc_arg == "1") {
    "oracle"
  } else {
    "aft_lognormal"
  }

} else {

  sc_methods <- trimws(
    unlist(
      strsplit(
        sc_arg,
        ",",
        fixed = TRUE
      )
    )
  )
}


valid_sc_methods <- c(
  "oracle",
  "rsf",
  "aft_lognormal",
  "km_x1",
  "km",
  "power_oracle"
)


if (any(!sc_methods %in% valid_sc_methods)) {

  stop(
    sprintf(
      "Unknown censoring model(s): %s",
      paste(
        setdiff(
          sc_methods,
          valid_sc_methods
        ),
        collapse = ", "
      )
    )
  )
}


only_cams <- if (length(args) >= 4L) {

  as.logical(
    as.integer(args[4])
  )

} else {

  NA
}


if (is.na(only_cams)) {
  only_cams <- FALSE
}


augmentation_arg <- if (length(args) >= 5L) {
  args[5]
} else {
  "wrong"
}


augmentation_methods <- trimws(
  unlist(
    strsplit(
      augmentation_arg,
      ",",
      fixed = TRUE
    )
  )
)


valid_augmentation_methods <- c(
  "same",
  "correct",
  "wrong",
  "oracle_event"
)


if (
  any(
    !augmentation_methods %in%
      valid_augmentation_methods
  )
) {

  stop(
    sprintf(
      "Unknown augmentation model(s): %s",
      paste(
        setdiff(
          augmentation_methods,
          valid_augmentation_methods
        ),
        collapse = ", "
      )
    )
  )
}


homoscedastic_event <- if (
  length(args) >= 6L
) {

  as.logical(
    as.integer(args[6])
  )

} else {

  FALSE
}


if (is.na(homoscedastic_event)) {
  homoscedastic_event <- FALSE
}


sc_ntree <- if (length(args) >= 7L) {
  as.integer(args[7])
} else {
  1000L
}


if (
  is.na(sc_ntree) ||
  sc_ntree < 1L
) {
  sc_ntree <- 1000L
}


num_runs <- if (length(args) >= 8L) {
  as.integer(args[8])
} else {
  10L
}


if (
  is.na(num_runs) ||
  num_runs < 1L
) {
  num_runs <- 10L
}


use_intersectional_R <- if (
  length(args) >= 9L
) {

  as.logical(
    as.integer(args[9])
  )

} else {

  FALSE
}


if (is.na(use_intersectional_R)) {
  use_intersectional_R <- FALSE
}


########################################
## Load libraries
########################################

suppressPackageStartupMessages(
  library(tidyverse)
)

suppressPackageStartupMessages(
  library(survival)
)

suppressPackageStartupMessages(
  library(quantreg)
)

suppressPackageStartupMessages(
  library(GauPro)
)

suppressPackageStartupMessages(
  library(gbm)
)

suppressPackageStartupMessages(
  library(grf)
)

suppressPackageStartupMessages(
  library(parallel)
)

suppressPackageStartupMessages(
  library(snow)
)


########################################
## Source code
########################################

source("./source_code.R")
source("./model_script.R")
source("./simu.R")
source("./fig.R")


########################################
## Global simulation configuration
########################################

alpha <- 0.1

xmin <- -2
xmax <- 2


gamma_list <- c(
  0.5,
  1,
  2,
  3
)


########################################
## Detect available cores
########################################

slurm_cores <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "SLURM_CPUS_PER_TASK"
    )
  )
)


num_cores <- if (
  is.na(slurm_cores) ||
  slurm_cores < 1L
) {

  detectCores()

} else {

  slurm_cores
}


cat(
  sprintf(
    paste0(
      "Starting sequential outer loop. ",
      "Inner algorithms will utilize ",
      "%d cores...\n"
    ),
    num_cores
  )
)


cat(
  "\nPer-setting configurations:\n"
)

print(
  setting_configs,
  row.names = FALSE
)


########################################
## Helpers
########################################

get_base_results_dir <- function(
  prob,
  only_cams
) {

  prob_label <- format(
    prob,
    trim = TRUE,
    scientific = FALSE
  )


  if (only_cams) {

    paste0(
      "../new_results",
      prob_label
    )

  } else {

    paste0(
      "../results",
      prob_label
    )
  }
}


event_scale_label <- if (
  homoscedastic_event
) {

  "event_constant_scale"

} else {

  "event_group_varying_scale"
}


########################################
## Run simulations
########################################

for (i in seq_len(num_runs)) {

  current_seed <- seed + (i - 1L)

  run_start_time <- proc.time()[3]


  for (sc_method in sc_methods) {

    for (
      augmentation_method
      in augmentation_methods
    ) {

      for (
        j in seq_len(
          nrow(setting_configs)
        )
      ) {

        current_config <-
          setting_configs[
            j,
            ,
            drop = FALSE
          ]


        setting <-
          current_config$setting[[1]]

        n_train <-
          current_config$n_train[[1]]

        n_calib <-
          current_config$n_calib[[1]]

        n_test <-
          current_config$n_test[[1]]

        bernoulli_prob <-
          current_config$bernoulli_prob[[1]]


        base_results_dir <-
          get_base_results_dir(
            bernoulli_prob,
            only_cams
          )


        run_folder0 <- file.path(
          base_results_dir,
          paste0(
            "sc_",
            sc_method
          ),
          paste0(
            "aug_",
            augmentation_method
          )
        )


        cat(
          sprintf(
            paste0(
              "\nConfiguration for %s:\n",
              "  n_train = %d\n",
              "  n_calib = %d\n",
              "  n_test = %d\n",
              "  bernoulli_prob = %.4f\n"
            ),
            setting,
            n_train,
            n_calib,
            n_test,
            bernoulli_prob
          )
        )


        ####################################
        ## Power-oracle censoring
        ####################################

        if (
          sc_method == "power_oracle"
        ) {

          for (gamma in gamma_list) {

            start_time <-
              proc.time()[3]


            cat(
              sprintf(
                paste0(
                  "\n=== Run %d | ",
                  "Setting: %s | ",
                  "G: %s | ",
                  "augmentation: %s | ",
                  "gamma: %.2f ===\n"
                ),
                i,
                setting,
                sc_method,
                augmentation_method,
                gamma
              )
            )


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

              use_oracle_sc =
                identical(
                  sc_method,
                  "oracle"
                ),

              sc_method =
                sc_method,

              augmentation_method =
                augmentation_method,

              homoscedastic_event =
                homoscedastic_event,

              sc_ntree =
                sc_ntree,

              gamma =
                gamma,

              use_intersectional_R =
                use_intersectional_R
            )


            run_folder_final <-
              file.path(
                run_folder0,
                paste0(
                  "gamma_",
                  gamma
                ),
                event_scale_label,
                as.character(
                  current_seed
                )
              )


            if (
              !dir.exists(
                run_folder_final
              )
            ) {

              dir.create(
                run_folder_final,
                showWarnings = FALSE,
                recursive = TRUE
              )
            }


            save_dir <- file.path(
              run_folder_final,
              sprintf(
                "%d. %s_seed_%d.csv",
                j,
                setting,
                current_seed
              )
            )


            write.csv(
              simures,
              save_dir,
              row.names = FALSE
            )


            cat(
              sprintf(
                paste0(
                  "%s for run %d ",
                  "(%s/%s/gamma=%.2f): ",
                  "%.2f seconds.\n"
                ),
                setting,
                i,
                sc_method,
                augmentation_method,
                gamma,
                proc.time()[3] -
                  start_time
              )
            )
          }


        ####################################
        ## Other censoring models
        ####################################

        } else {

          start_time <-
            proc.time()[3]


          cat(
            sprintf(
              paste0(
                "\n=== Run %d | ",
                "Setting: %s | ",
                "G: %s | ",
                "augmentation: %s ===\n"
              ),
              i,
              setting,
              sc_method,
              augmentation_method
            )
          )


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

            use_oracle_sc =
              identical(
                sc_method,
                "oracle"
              ),

            sc_method =
              sc_method,

            augmentation_method =
              augmentation_method,

            homoscedastic_event =
              homoscedastic_event,

            sc_ntree =
              sc_ntree,

            use_intersectional_R =
              use_intersectional_R
          )


          run_folder_final <-
            file.path(
              run_folder0,
              event_scale_label,
              as.character(
                current_seed
              )
            )


          if (
            !dir.exists(
              run_folder_final
            )
          ) {

            dir.create(
              run_folder_final,
              showWarnings = FALSE,
              recursive = TRUE
            )
          }


          save_dir <- file.path(
            run_folder_final,
            sprintf(
              "%d. %s_seed_%d.csv",
              j,
              setting,
              current_seed
            )
          )


          write.csv(
            simures,
            save_dir,
            row.names = FALSE
          )


          cat(
            sprintf(
              paste0(
                "%s for run %d ",
                "(%s/%s): ",
                "%.2f seconds.\n"
              ),
              setting,
              i,
              sc_method,
              augmentation_method,
              proc.time()[3] -
                start_time
            )
          )
        }
      }
    }
  }


  cat(
    sprintf(
      paste0(
        "\nCompleted run %d/%d ",
        "in %.2f seconds.\n"
      ),
      i,
      num_runs,
      proc.time()[3] -
        run_start_time
    )
  )
}


########################################
## Create plots and summaries
########################################

bernoulli_probs <- unique(
  setting_configs$bernoulli_prob
)


for (
  bernoulli_prob
  in bernoulli_probs
) {

  base_results_dir <-
    get_base_results_dir(
      bernoulli_prob,
      only_cams
    )


  plots_root_dir <- sub(
    "results",
    "plots",
    base_results_dir,
    fixed = TRUE
  )


  summaries_root_dir <- sub(
    "results",
    "summaries",
    base_results_dir,
    fixed = TRUE
  )


  for (
    sc_method
    in sc_methods
  ) {

    for (
      augmentation_method
      in augmentation_methods
    ) {

      method_results_dir <- file.path(
        base_results_dir,
        paste0(
          "sc_",
          sc_method
        ),
        paste0(
          "aug_",
          augmentation_method
        )
      )


      method_plots_dir <- file.path(
        plots_root_dir,
        paste0(
          "sc_",
          sc_method
        ),
        paste0(
          "aug_",
          augmentation_method
        )
      )


      method_summaries_dir <- file.path(
        summaries_root_dir,
        paste0(
          "sc_",
          sc_method
        ),
        paste0(
          "aug_",
          augmentation_method
        )
      )


      ####################################
      ## Power oracle
      ####################################

      if (
        sc_method == "power_oracle"
      ) {

        for (
          gamma
          in gamma_list
        ) {

          gamma_label <- paste0(
            "gamma_",
            format(
              gamma,
              trim = TRUE,
              scientific = FALSE
            )
          )


          current_results_dir <-
            file.path(
              method_results_dir,
              gamma_label,
              event_scale_label
            )


          current_plots_dir <-
            file.path(
              method_plots_dir,
              gamma_label,
              event_scale_label
            )


          current_summaries_dir <-
            file.path(
              method_summaries_dir,
              gamma_label,
              event_scale_label
            )


          if (
            dir.exists(
              current_results_dir
            )
          ) {

            make_plots(
              results_dir =
                current_results_dir,

              plots_dir =
                current_plots_dir,

              summaries_dir =
                current_summaries_dir,

              target_alpha =
                alpha,

              use_intersectional_R =
                use_intersectional_R
            )

          } else {

            warning(
              sprintf(
                paste0(
                  "Skipping plots; ",
                  "results directory ",
                  "does not exist: %s"
                ),
                current_results_dir
              )
            )
          }
        }


      ####################################
      ## Other censoring models
      ####################################

      } else {

        current_results_dir <-
          file.path(
            method_results_dir,
            event_scale_label
          )


        current_plots_dir <-
          file.path(
            method_plots_dir,
            event_scale_label
          )


        current_summaries_dir <-
          file.path(
            method_summaries_dir,
            event_scale_label
          )


        if (
          dir.exists(
            current_results_dir
          )
        ) {

          make_plots(
            results_dir =
              current_results_dir,

            plots_dir =
              current_plots_dir,

            summaries_dir =
              current_summaries_dir,

            target_alpha =
              alpha,

            use_intersectional_R =
              use_intersectional_R
          )

        } else {

          warning(
            sprintf(
              paste0(
                "Skipping plots; ",
                "results directory ",
                "does not exist: %s"
              ),
              current_results_dir
            )
          )
        }
      }
    }
  }
}


cat(
  sprintf(
    "\nCompleted in %.2f seconds.\n",
    proc.time()[3] -
      total_start_time
  )
)