#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
repo_root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
lib <- Sys.getenv("R_LIBS_USER", file.path(path.expand("~"), "R/support2-library"))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
needed <- c("survival", "quantreg")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = lib, repos = "https://cloud.r-project.org")
if (!all(vapply(needed, requireNamespace, logical(1), quietly = TRUE))) stop("Package setup failed.")
source(file.path(repo_root, "real_data/support2.R"))
raw <- load_support2(cache_dir = file.path(repo_root, "real_data/data"))
cat("Ready:", nrow(raw), "patients; CSV MD5:", attr(raw, "source_md5"), "\n")
cat("R packages:", paste(needed, vapply(needed, function(p) as.character(packageVersion(p)), "")), "\n")
