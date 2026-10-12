# Run from the project root. This permanently deletes native folders and their contents.
tasks_dir <- file.path("results", "support2", "experiment_A", "tasks")
if (!dir.exists(tasks_dir)) {
  stop("Tasks directory not found: ", tasks_dir)
}

seed_dirs <- list.dirs(tasks_dir, full.names = TRUE, recursive = FALSE)
seed_dirs <- seed_dirs[grepl("^seed_[0-9]+$", basename(seed_dirs))]
native_dirs <- file.path(seed_dirs, "native")
native_dirs <- native_dirs[dir.exists(native_dirs)]

for (native_dir in native_dirs) {
  unlink(native_dir, recursive = TRUE)
  if (dir.exists(native_dir)) {
    stop("Failed to delete: ", native_dir)
  }
  message("Deleted: ", native_dir)
}

message("Deleted ", length(native_dirs), " native folder(s).")
