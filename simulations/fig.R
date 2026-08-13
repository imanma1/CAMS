make_plots <- function(
    results_dir = "../results",
    plots_dir = "../plots",
    summaries_dir = "../summaries",
    target_alpha = 0.1,
    use_intersectional_R = FALSE
) {

  target_cov <- 1 - target_alpha

  dir.create(
    plots_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  dir.create(
    summaries_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  # ---------------------------------------------------------
  # 1. Detect seed folders
  # ---------------------------------------------------------
  seed_dirs_full <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
  seed_dir_names <- basename(seed_dirs_full)

  is_seed_dir <- grepl("^[0-9]+$", seed_dir_names)

  seed_dirs_full <- seed_dirs_full[is_seed_dir]
  seed_dir_names <- seed_dir_names[is_seed_dir]

  if (length(seed_dirs_full) == 0) {
    stop(sprintf("No numeric seed folders found in %s", results_dir))
  }

  seed_order <- order(as.integer(seed_dir_names))
  seed_dirs_full <- seed_dirs_full[seed_order]
  seed_dir_names <- seed_dir_names[seed_order]

  cat(sprintf("Detected %d seed folders.\n", length(seed_dirs_full)))

  # ---------------------------------------------------------
  # 2. Detect settings from first seed folder
  #    Expected format:
  #    "1. homo_cens_seed_1.csv"
  # ---------------------------------------------------------
  first_seed_dir <- seed_dirs_full[1]

  csv_files_first_seed <- list.files(
    first_seed_dir,
    pattern = "\\.csv$",
    full.names = FALSE
  )

  if (length(csv_files_first_seed) == 0) {
    stop(sprintf("No CSV files found in first seed folder: %s", first_seed_dir))
  }

  file_info <- data.frame(
    file = csv_files_first_seed,
    stringsAsFactors = FALSE
  )

  file_info$order <- as.integer(sub("^([0-9]+)\\.\\s+.*$", "\\1", file_info$file))

  file_info$setting <- sub("^[0-9]+\\.\\s+", "", file_info$file)
  file_info$setting <- sub("_seed_[0-9]+\\.csv$", "", file_info$setting)

  file_info <- file_info[!is.na(file_info$order), , drop = FALSE]

  if (nrow(file_info) == 0) {
    stop("No CSV files matched the expected format: '<order>. <setting>_seed_<seed>.csv'")
  }

  file_info <- file_info[order(file_info$order), , drop = FALSE]

  cat("Detected settings in this order:\n")
  print(file_info[, c("order", "setting")])

  # ---------------------------------------------------------
  # 3. Helper for safe y-limits
  # ---------------------------------------------------------
  get_ylim <- function(values) {
    y_min <- suppressWarnings(min(values, na.rm = TRUE))
    y_max <- suppressWarnings(max(values, na.rm = TRUE))

    if (is.infinite(y_min) || is.infinite(y_max) || is.na(y_min) || is.na(y_max)) {
      return(c(0, 1))
    }

    if (y_min == y_max) {
      padding <- ifelse(abs(y_min) < 1e-8, 0.1, abs(y_min) * 0.1)
      return(c(y_min - padding, y_max + padding))
    }

    padding <- 0.05 * (y_max - y_min)
    c(y_min - padding, y_max + padding)
  }

  xml_escape <- function(x) {

    x <- as.character(x)

    x[is.na(x)] <- ""

    x <- gsub(
      "&",
      "&amp;",
      x,
      fixed = TRUE
    )

    x <- gsub(
      "<",
      "&lt;",
      x,
      fixed = TRUE
    )

    x <- gsub(
      ">",
      "&gt;",
      x,
      fixed = TRUE
    )

    x <- gsub(
      "\"",
      "&quot;",
      x,
      fixed = TRUE
    )

    x <- gsub(
      "'",
      "&apos;",
      x,
      fixed = TRUE
    )

    x
  }


  write_summary_xml <- function(
      summary_df,
      file_path,
      setting
  ) {

    xml_lines <- c(
      '<?xml version="1.0" encoding="UTF-8"?>',
      sprintf(
        '<simulation_summary setting="%s">',
        xml_escape(setting)
      )
    )

    for (row_idx in seq_len(nrow(summary_df))) {

      method_name <- summary_df$method[
        row_idx
      ]

      xml_lines <- c(
        xml_lines,
        sprintf(
          '  <method name="%s">',
          xml_escape(method_name)
        )
      )

      metric_names <- setdiff(
        colnames(summary_df),
        "method"
      )

      for (metric_name in metric_names) {

        value <- summary_df[
          row_idx,
          metric_name
        ][[1]]

        if (
          length(value) == 0L ||
            is.na(value)
        ) {
          value <- ""
        } else {
          value <- as.character(value)
        }

        xml_lines <- c(
          xml_lines,
          sprintf(
            '    <metric name="%s">%s</metric>',
            xml_escape(metric_name),
            xml_escape(value)
          )
        )
      }

      xml_lines <- c(
        xml_lines,
        "  </method>"
      )
    }

    xml_lines <- c(
      xml_lines,
      "</simulation_summary>"
    )

    writeLines(
      xml_lines,
      con = file_path,
      useBytes = TRUE
    )
  }

  # ---------------------------------------------------------
  # 4. Plotting helper
  # ---------------------------------------------------------
  draw_paper_boxplot <- function(data,
                                 metric_col,
                                 y_label,
                                 add_target_line = FALSE,
                                 method_levels) {

    ylim <- get_ylim(data[[metric_col]])

    n_methods <- length(method_levels)
    border_cols <- rep(seq_len(max(8, n_methods)) + 1, length.out = n_methods)

    boxplot(
      data[[metric_col]] ~ data$method,
      col = "white",
      border = border_cols,
      ylab = y_label,
      xlab = "",
      las = 2,
      cex.axis = 0.85,
      ylim = ylim,
      main = ""
    )

    if (add_target_line) {
      abline(h = target_cov, lty = 2, col = "black", lwd = 1.5)
    }
  }

  # ---------------------------------------------------------
  # 5. Metric columns
  # ---------------------------------------------------------
  if (use_intersectional_R) {

    metric_specs <- list(

      list(
        col = "Marginal coverage",
        label = "Marginal coverage",
        target = TRUE
      ),

      list(
        col = paste0(
          "group coverage for ",
          "x_1 = 0 and x_2 <= 0"
        ),
        label = "Coverage: X1 = 0, X2 <= 0",
        target = TRUE
      ),

      list(
        col = paste0(
          "group coverage for ",
          "x_1 = 0 and x_2 > 0"
        ),
        label = "Coverage: X1 = 0, X2 > 0",
        target = TRUE
      ),

      list(
        col = paste0(
          "group coverage for ",
          "x_1 = 1 and x_2 <= 0"
        ),
        label = "Coverage: X1 = 1, X2 <= 0",
        target = TRUE
      ),

      list(
        col = paste0(
          "group coverage for ",
          "x_1 = 1 and x_2 > 0"
        ),
        label = "Coverage: X1 = 1, X2 > 0",
        target = TRUE
      ),

      list(
        col = "lower bound values mean",
        label = "Overall average lower bound",
        target = FALSE
      ),

      list(
        col = paste0(
          "lower bound mean for ",
          "x_1 = 0 and x_2 <= 0"
        ),
        label = "Lower bound: X1 = 0, X2 <= 0",
        target = FALSE
      ),

      list(
        col = paste0(
          "lower bound mean for ",
          "x_1 = 0 and x_2 > 0"
        ),
        label = "Lower bound: X1 = 0, X2 > 0",
        target = FALSE
      ),

      list(
        col = paste0(
          "lower bound mean for ",
          "x_1 = 1 and x_2 <= 0"
        ),
        label = "Lower bound: X1 = 1, X2 <= 0",
        target = FALSE
      ),

      list(
        col = paste0(
          "lower bound mean for ",
          "x_1 = 1 and x_2 > 0"
        ),
        label = "Lower bound: X1 = 1, X2 > 0",
        target = FALSE
      )
    )

  } else {

    metric_specs <- list(

      list(
        col = "Marginal coverage",
        label = "Marginal coverage",
        target = TRUE
      ),

      list(
        col = "group coverage for x_1 = 0",
        label = "Coverage: X1 = 0",
        target = TRUE
      ),

      list(
        col = "group coverage for x_1 = 1",
        label = "Coverage: X1 = 1",
        target = TRUE
      ),

      list(
        col = "lower bound values mean",
        label = "Overall average lower bound",
        target = FALSE
      ),

      list(
        col = "lower bound mean for x_1 = 0",
        label = "Lower bound: X1 = 0",
        target = FALSE
      ),

      list(
        col = "lower bound mean for x_1 = 1",
        label = "Lower bound: X1 = 1",
        target = FALSE
      )
    )
  }

  # ---------------------------------------------------------
  # 6. Loop through detected settings
  # ---------------------------------------------------------
  for (row_id in seq_len(nrow(file_info))) {

    plot_order <- file_info$order[row_id]
    setting <- file_info$setting[row_id]

    cat(sprintf("\nProcessing setting: %s...\n", setting))

    all_data <- data.frame()
    missing_files_count <- 0

    for (k in seq_along(seed_dirs_full)) {

      seed_dir <- seed_dirs_full[k]
      seed_name <- seed_dir_names[k]

      # New filename format:
      # "<order>. <setting>_seed_<seed>.csv"
      file_name <- sprintf("%d. %s_seed_%s.csv",
                           plot_order, setting, seed_name)

      file_path <- file.path(seed_dir, file_name)

      if (file.exists(file_path)) {

        temp_df <- read.csv(
          file_path,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )

        temp_df$source_seed <- seed_name
        all_data <- rbind(all_data, temp_df)

      } else {

        # Fallback: search by setting + seed in case the prefix differs
        pattern <- sprintf("^[0-9]+\\.\\s+%s_seed_%s\\.csv$",
                           setting, seed_name)

        candidate_files <- list.files(
          seed_dir,
          pattern = pattern,
          full.names = TRUE
        )

        if (length(candidate_files) == 1) {

          temp_df <- read.csv(
            candidate_files[1],
            check.names = FALSE,
            stringsAsFactors = FALSE
          )

          temp_df$source_seed <- seed_name
          all_data <- rbind(all_data, temp_df)

        } else {
          cat(sprintf("  -> [WARNING] Missing file skipped: %s\n", file_path))
          missing_files_count <- missing_files_count + 1
        }
      }
    }

    if (nrow(all_data) == 0) {
      cat(sprintf("  -> [ERROR] No data found for %s. Skipping.\n", setting))
      next
    }

    if (!("method" %in% colnames(all_data))) {
      cat(sprintf("  -> [ERROR] Column 'method' not found for %s. Skipping.\n", setting))
      next
    }

    missing_metric_cols <- vapply(
      metric_specs,
      function(spec) !(spec$col %in% colnames(all_data)),
      logical(1)
    )

    if (any(missing_metric_cols)) {
      missing_names <- vapply(
        metric_specs[missing_metric_cols],
        function(spec) spec$col,
        character(1)
      )

      cat(sprintf(
        "  -> [ERROR] Missing metric columns for %s: %s. Skipping.\n",
        setting,
        paste(missing_names, collapse = ", ")
      ))

      next
    }

    if (missing_files_count > 0) {
      cat(sprintf(
        "  -> Plotting %s using %d available files out of %d seed folders.\n",
        setting,
        length(seed_dirs_full) - missing_files_count,
        length(seed_dirs_full)
      ))
    }

    # Method order is detected from the CSV rows
    method_levels <- unique(all_data$method)
    all_data$method <- factor(all_data$method, levels = method_levels)

    cat("  -> Detected methods in this order:\n")
    print(method_levels)

    # ---------------------------------------------------------
    # Report simulation means while explicitly tracking failures.
    #
    # A method-run is considered complete only if both its
    # marginal coverage and overall lower-bound mean are finite.
    #
    # This prevents numerical/selection failures from being
    # silently hidden by na.rm = TRUE.
    # ---------------------------------------------------------

    numeric_cols <- names(all_data)[
      vapply(
        all_data,
        is.numeric,
        logical(1)
      )
    ]

    summary_rows <- lapply(
      method_levels,
      function(method_name) {

        method_data <- all_data[
          all_data$method == method_name,
          ,
          drop = FALSE
        ]

        complete_run <- (
          is.finite(
            method_data[["Marginal coverage"]]
          ) &
          is.finite(
            method_data[["lower bound values mean"]]
          )
        )

        available_seeds <- unique(
          method_data$source_seed
        )

        complete_seeds <- unique(
          method_data$source_seed[
            complete_run
          ]
        )

        failed_seeds <- unique(
          method_data$source_seed[
            !complete_run
          ]
        )

        n_available <- length(
          available_seeds
        )

        n_complete <- length(
          complete_seeds
        )

        n_failed <- length(
          failed_seeds
        )

        failure_rate <- if (
          n_available > 0L
        ) {
          n_failed / n_available
        } else {
          NA_real_
        }

        row <- data.frame(
          method = method_name,
          `number of available seeds` =
            n_available,
          `number of complete seeds` =
            n_complete,
          `number of failed seeds` =
            n_failed,
          `failure rate` =
            failure_rate,
          `failed seed IDs` =
            if (n_failed > 0L) {
              paste(
                failed_seeds,
                collapse = ","
              )
            } else {
              ""
            },
          check.names = FALSE,
          stringsAsFactors = FALSE
        )

        for (metric_name in numeric_cols) {

          values <- method_data[[metric_name]]

          complete_values <- values[complete_run]

          finite_values <- complete_values[
            is.finite(complete_values)
          ]

          row[[
            paste0(
              metric_name,
              " mean among complete seeds"
            )
          ]] <- if (length(finite_values) > 0L) {
            mean(finite_values)
          } else {
            NA_real_
          }

          row[[
            paste0(
              metric_name,
              " sd across complete seeds"
            )
          ]] <- if (length(finite_values) >= 2L) {
            stats::sd(finite_values)
          } else {
            NA_real_
          }
        }

        row
      }
    )

    summary_df <- do.call(
      rbind,
      summary_rows
    )
    summary_csv_name <- file.path(
      summaries_dir,
      sprintf(
        "%d. %s_summary.csv",
        plot_order,
        setting
      )
    )

    summary_xml_name <- file.path(
      summaries_dir,
      sprintf(
        "%d. %s_summary.xml",
        plot_order,
        setting
      )
    )

    write.csv(
      summary_df,
      summary_csv_name,
      row.names = FALSE
    )

    write_summary_xml(
      summary_df = summary_df,
      file_path = summary_xml_name,
      setting = setting
    )

    cat(
      sprintf(
        "  -> CSV summary saved to: %s\n",
        summary_csv_name
      )
    )

    cat(
      sprintf(
        "  -> XML summary saved to: %s\n",
        summary_xml_name
      )
    )

    # Output filename keeps your desired format:
    # "<order>. <setting>.png"
    png_name <- file.path(plots_dir, sprintf("%d. %s.png", plot_order, setting))

    if (use_intersectional_R) {

      plot_rows <- 3
      plot_columns <- 4
      plot_width <- 20
      plot_height <- 14

    } else {

      plot_rows <- 2
      plot_columns <- 3
      plot_width <- 15
      plot_height <- 11
    }

    png(
      filename = png_name,
      width = plot_width,
      height = plot_height,
      units = "in",
      res = 300
    )

    par(
      mfrow = c(
        plot_rows,
        plot_columns
      ),
      mar = c(11, 4, 2, 1),
      oma = c(0, 0, 3, 0)
    )

    for (spec in metric_specs) {
      draw_paper_boxplot(
        data = all_data,
        metric_col = spec$col,
        y_label = spec$label,
        add_target_line = spec$target,
        method_levels = method_levels
      )
    }

    number_of_empty_panels <- (
      plot_rows * plot_columns
    ) - length(metric_specs)

    if (number_of_empty_panels > 0L) {
      for (
        empty_panel in seq_len(
          number_of_empty_panels
        )
      ) {
        plot.new()
      }
    }

    mtext(
      sprintf("Results for Setting: %s", setting),
      side = 3,
      outer = TRUE,
      cex = 1.8,
      font = 2
    )

    dev.off()

    cat(sprintf("  -> Plot saved to: %s\n", png_name))
  }

  cat("\nAll detected settings processed successfully!\n")
}
