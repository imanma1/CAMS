#!/usr/bin/env Rscript
# Merge isolated native Experiment A tasks after the array finishes.
args <- commandArgs(trailingOnly = TRUE)
opts <- list(input = "", output = "", `expected-seeds` = "")
if ("--help" %in% args) {
  cat("Rscript real_data/merge_support2.R --input=/path/to/experiment_A_JOBID --expected-seeds=1:100\n")
  quit(status = 0)
}
for (arg in args) {
  pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  if (!grepl("^--[^=]+=.+$", arg) || !pair[1] %in% names(opts)) stop("Invalid argument: ", arg)
  opts[[pair[1]]] <- paste(pair[-1], collapse = "=")
}
if (!nzchar(opts$input) || !dir.exists(opts$input)) stop("Provide an existing --input directory.")
file_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
repo <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
source(file.path(repo, "real_data/support2.R"))
task_parent <- file.path(opts$input, "tasks")
tasks <- if (dir.exists(task_parent)) sort(list.dirs(task_parent, recursive = FALSE)) else opts$input
if (!length(tasks)) stop("No task directories found.")
metrics <- list(); statuses <- list(); method_status <- list(); signature <- NULL
for (task in tasks) {
  protocol_path <- file.path(task, "study_config.R")
  if (file.exists(protocol_path)) {
    protocol <- dget(protocol_path)
    if (is.null(signature)) signature <- protocol else if (!identical(signature, protocol))
      stop("Task settings/code/data differ; refusing to combine: ", task)
  }
  status_path <- file.path(task, "run_status.csv")
  if (!file.exists(status_path)) {
    seed <- suppressWarnings(as.integer(sub("^seed_", "", basename(task))))
    statuses[[length(statuses) + 1L]] <- data.frame(seed = seed, scenario = "native",
      status = "incomplete", elapsed_seconds = NA_real_, message = "No finalized run status")
    next
  }
  status <- read.csv(status_path, check.names = FALSE)
  if (any(status$scenario != "native")) stop("This merge command is for native Experiment A only.")
  for (ii in seq_len(nrow(status))) {
    run_dir <- file.path(task, "native", paste0("seed_", status$seed[ii]))
    if (status$status[ii] %in% c("ok", "partial")) {
      required <- c("COMPLETE", "config.R", "metrics.csv", "method_status.csv")
      if (!all(file.exists(file.path(run_dir, required)))) {
        status$status[ii] <- "incomplete"; status$message[ii] <- "Missing finalized per-run output"
        next
      }
      m <- read.csv(file.path(run_dir, "metrics.csv"), check.names = FALSE)
      if (any(m$seed != status$seed[ii]) || any(m$scenario != "native")) stop("Mismatched run labels.")
      metrics[[length(metrics) + 1L]] <- m
      ms <- read.csv(file.path(run_dir, "method_status.csv"), check.names = FALSE)
      ms$seed <- status$seed[ii]; ms$scenario <- "native"
      method_status[[length(method_status) + 1L]] <- ms
      if (any(ms$status != "ok")) status$status[ii] <- "partial"
    }
  }
  statuses[[length(statuses) + 1L]] <- status
}
runs <- do.call(rbind, statuses)
if (anyDuplicated(paste(runs$scenario, runs$seed))) stop("Duplicate seeds; refusing to double-count splits.")
expected <- integer()
if (nzchar(opts[["expected-seeds"]])) {
  text <- opts[["expected-seeds"]]
  if (grepl("^[0-9]+:[0-9]+$", text)) {
    edges <- as.integer(strsplit(text, ":", fixed = TRUE)[[1]])
    if (anyNA(edges) || edges[2] < edges[1]) stop("Invalid expected seed range.")
    expected <- seq.int(edges[1], edges[2])
  } else {
    tokens <- strsplit(text, ",", fixed = TRUE)[[1]]
    if (any(!grepl("^[0-9]+$", tokens))) stop("Invalid expected seed list.")
    expected <- as.integer(tokens)
    if (anyNA(expected) || anyDuplicated(expected)) stop("Invalid expected seed list.")
  }
  extra <- setdiff(runs$seed, expected)
  if (length(extra)) stop("Unexpected seeds: ", paste(extra, collapse = ","))
  missing <- setdiff(expected, runs$seed)
  if (length(missing)) runs <- rbind(runs, data.frame(seed = missing, scenario = "native",
    status = "missing", elapsed_seconds = NA_real_, message = "Expected array task has no results"))
}
out <- if (nzchar(opts$output)) opts$output else file.path(opts$input, "merged")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(runs[order(runs$seed), ], file.path(out, "run_status.csv"), row.names = FALSE)
if (length(metrics)) {
  m <- do.call(rbind, metrics)
  keys <- paste(m$scenario, m$seed, m$method, m$audit_group, sep = "\r")
  if (anyDuplicated(keys)) stop("Duplicate metric rows.")
  expected_methods <- unique(m$method)
  for (s in unique(m$seed)) {
    mm <- m[m$seed == s, ]
    if (!setequal(mm$method, expected_methods)) stop("Inconsistent method lists across seeds.")
  }
  write.csv(m, file.path(out, "metrics_all.csv"), row.names = FALSE)
  write.csv(summarize_support2(m), file.path(out, "summary.csv"), row.names = FALSE)
  finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  finite_sd <- function(x) if (sum(is.finite(x)) > 1) sd(x[is.finite(x)]) else NA_real_
  cells <- unique(m[, c("scenario", "method", "audit_group")])
  by_group <- lapply(seq_len(nrow(cells)), function(ii) {
    x <- m[m$scenario == cells$scenario[ii] & m$method == cells$method[ii] &
      m$audit_group == cells$audit_group[ii], ]
    cbind(cells[ii, ], data.frame(n_runs = nrow(x), n_complete = sum(x$complete),
      n_exact = sum(is.finite(x$coverage_exact)),
      mean_exact_coverage = finite_mean(x$coverage_exact), sd_exact_coverage = finite_sd(x$coverage_exact),
      mean_ipcw_coverage = finite_mean(x$coverage_ipcw), mean_aipcw_coverage = finite_mean(x$coverage_aipcw),
      mean_bound = finite_mean(x$mean_bound), sd_bound = finite_sd(x$mean_bound)))
  })
  write.csv(do.call(rbind, by_group), file.path(out, "summary_by_group.csv"), row.names = FALSE)
}
if (length(method_status)) write.csv(do.call(rbind, method_status),
  file.path(out, "method_status_all.csv"), row.names = FALSE)
if (!is.null(signature)) dput(signature, file.path(out, "study_config.R"))
writeLines(c("Repeated splits share patients; split standard deviations are descriptive, not independent-cohort standard errors.",
  "DFT uses conditional imputation of missing censoring times, not known native Type-I times.",
  "Exact coverage is unavailable when any original coverage indicator in that audit group is unresolved.",
  "Missing/failed tasks are reported in run_status.csv; completion counts accompany means."), file.path(out, "INTERPRETATION.txt"))
cat("Merged results:", normalizePath(out), "\n")
if (any(runs$status != "ok")) {
  message("Missing, incomplete, or failed method comparisons remain; inspect run_status.csv.")
  quit(status = 1)
}
