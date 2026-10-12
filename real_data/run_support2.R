#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
opts <- list(scenario = "native", seeds = "1", horizon = "180", alpha = "0.1",
  `feature-set` = "entry", `landmark-days` = "0", `probability-floor` = "0.05",
  `grid-step` = "0.001", `integration-steps` = "256", `include-subgroups` = "true",
  `include-cqr` = "true", `include-dft` = "true", `dft-fixed-cutoff` = "auto",
  `imputation-seed` = "", `strict-methods` = "true", output = "", csv = "")
opts[["aft-ridge-fallback"]] <- "true"
opts[["aft-ridge-theta"]] <- "0.1"
if ("--help" %in% args) {
  cat("Run from any directory: Rscript /path/to/CAMS/real_data/run_support2.R [options]\n",
    "--scenario=all|native|added_common|added_group --seeds=1,2,3 --horizon=180\n",
    "--alpha=0.1 --feature-set=entry|day3 --landmark-days=0 --probability-floor=0.05\n",
    "--grid-step=0.001 --integration-steps=256 --include-subgroups=true --include-cqr=true\n",
    "--include-dft=true --dft-fixed-cutoff=auto|days --imputation-seed=integer\n",
    "--strict-methods=true (save partial results and exit nonzero on a method failure)\n",
    "--aft-ridge-fallback=true --aft-ridge-theta=0.1\n",
    "--csv=/path/to/original/support2.csv --output=/path/to/new/results\n", sep = "")
  quit(status = 0)
}
for (arg in args) {
  if (!grepl("^--[^=]+=.+$", arg)) stop("Use --name=value arguments: ", arg)
  pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  name <- pair[1]; value <- paste(pair[-1], collapse = "=")
  if (!name %in% names(opts)) stop("Unknown option: ", name)
  opts[[name]] <- value
}
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
repo_root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
source(file.path(repo_root, "real_data/support2.R"))
as_bool <- function(x) {
  if (!x %in% c("true", "false")) stop("Boolean options must be true or false.")
  identical(x, "true")
}
seeds <- as.integer(strsplit(opts$seeds, ",", fixed = TRUE)[[1]])
seed_tokens <- strsplit(opts$seeds, ",", fixed = TRUE)[[1]]
if (!length(seeds) || any(!grepl("^[0-9]+$", seed_tokens)) || anyNA(seeds) ||
    any(seeds < 0) || anyDuplicated(seeds)) stop("Invalid seeds.")
strict_methods <- as_bool(opts[["strict-methods"]])
scenarios <- if (opts$scenario == "all") c("native", "added_common", "added_group") else opts$scenario
if (!all(scenarios %in% c("native", "added_common", "added_group"))) stop("Invalid scenario.")
output_root <- if (nzchar(opts$output)) opts$output else file.path(repo_root,
  "results/support2", paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
signature <- opts[setdiff(names(opts), c("seeds", "output", "csv", "strict-methods", "scenario"))]
signature$code_md5 <- unname(tools::md5sum(file.path(repo_root, c("utils/cams.R",
  "utils/cqr.R", "utils/utils.R", "real_data/support2.R", "real_data/baselines.R"))))
signature$source_md5 <- attr(load_support2(if (nzchar(opts$csv)) opts$csv else NULL,
  file.path(repo_root, "real_data/data")), "source_md5")
signature_path <- file.path(output_root, "study_config.R")
if (file.exists(signature_path)) {
  if (!identical(dget(signature_path), signature)) stop("Output root has different settings/code/data. Use a new output root.")
} else {
  if (file.exists(file.path(output_root, "metrics_all.csv"))) stop("Legacy output root: use a new directory for the expanded comparison.")
  dput(signature, signature_path)
}
all_metrics <- list(); runs <- list()
# Resume with new seeds without discarding previous run summaries.
metrics_path <- file.path(output_root, "metrics_all.csv")
status_path <- file.path(output_root, "run_status.csv")
if (file.exists(metrics_path)) all_metrics[[1L]] <- read.csv(metrics_path, check.names = FALSE)
if (file.exists(status_path)) runs[[1L]] <- read.csv(status_path, check.names = FALSE)
for (seed in seeds) for (scenario in scenarios) {
  destination <- file.path(output_root, scenario, paste0("seed_", seed))
  if (dir.exists(destination)) stop("Refusing to overwrite an existing run: ", destination)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  started <- proc.time()[3]
  run <- tryCatch(run_support2(repo_root = repo_root,
    csv_path = if (nzchar(opts$csv)) opts$csv else NULL, seed = seed, scenario = scenario,
    horizon = as.numeric(opts$horizon), alpha = as.numeric(opts$alpha),
    feature_set = opts[["feature-set"]], landmark_days = as.numeric(opts[["landmark-days"]]),
    probability_floor = as.numeric(opts[["probability-floor"]]),
    grid_step = as.numeric(opts[["grid-step"]]),
    integration_steps = as.integer(opts[["integration-steps"]]),
    include_subgroups = as_bool(opts[["include-subgroups"]]),
    include_cqr = as_bool(opts[["include-cqr"]]), output_dir = destination,
    include_dft = as_bool(opts[["include-dft"]]), dft_fixed_cutoff = opts[["dft-fixed-cutoff"]],
    imputation_seed = if (nzchar(opts[["imputation-seed"]]))
      as.numeric(opts[["imputation-seed"]]) else NULL,
    aft_ridge_fallback = as_bool(opts[["aft-ridge-fallback"]]),
    aft_ridge_theta = as.numeric(opts[["aft-ridge-theta"]])), error = function(e) e)
  failed <- inherits(run, "error")
  partial <- !failed && any(run$status$status != "ok")
  runs[[length(runs) + 1L]] <- data.frame(seed = seed, scenario = scenario,
    status = if (failed) "failed" else if (partial) "partial" else "ok",
    elapsed_seconds = unname(proc.time()[3] - started),
    message = if (failed) conditionMessage(run) else if (partial)
      paste(run$status$method[run$status$status != "ok"], collapse = "; ") else "")
  if (failed) {
    message("Failed: ", scenario, "/", seed, ": ", conditionMessage(run))
    writeLines(conditionMessage(run), file.path(destination, "error.txt"))
  } else {
    all_metrics[[length(all_metrics) + 1L]] <- run$metrics
    writeLines(if (partial) "partial" else "ok", file.path(destination, "COMPLETE"))
    if (partial) message("Partial comparison: inspect ", file.path(destination, "method_status.csv"))
  }
  write.csv(do.call(rbind, runs), file.path(output_root, "run_status.csv"), row.names = FALSE)
  if (length(all_metrics)) write.csv(do.call(rbind, all_metrics),
    file.path(output_root, "metrics_all.csv"), row.names = FALSE)
}
cat("Results: ", normalizePath(output_root), "\n", sep = "")
if (length(all_metrics)) write.csv(summarize_support2(do.call(rbind, all_metrics)),
  file.path(output_root, "summary.csv"), row.names = FALSE)
if (any(do.call(rbind, runs)$status == "failed") ||
    (strict_methods && any(do.call(rbind, runs)$status == "partial"))) quit(status = 1)
