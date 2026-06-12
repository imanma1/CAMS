make_plots <- function() {
  # 1. Define Parameters
  # ---------------------------------------------------------
  results_dir <- "../old_results"
  plots_dir <- file.path("..", "plots")
  num_runs <- 50             # The number of seed folders to loop through
  target_alpha <- 0.1          # Miscoverage target
  target_cov <- 1 - target_alpha

  # List of all settings used in the simulations
  setting_list <- c("homo_cens", "cov_cens", "prot_cens", 
                    "heavy_prot_cens", "heavy_inter_cens", 
                    "surv_misspec", "cens_misspec", "simul_misspec")

  # Create the plots directory if it doesn't exist
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  # 2. Define Plotting Helper Function
  # ---------------------------------------------------------
  # This function replicates the visual style of the paper's figures
  draw_paper_boxplot <- function(data, metric_col, y_label, add_target_line = FALSE) {
    
    # Calculate y-limits dynamically to add a little padding
    # (Fallback to 0-1 if data is empty or invalid to prevent plot crashes)
    y_min <- suppressWarnings(min(data[[metric_col]], na.rm = TRUE))
    y_max <- suppressWarnings(max(data[[metric_col]], na.rm = TRUE))
    if(is.infinite(y_min)) { y_min <- 0; y_max <- 1 }
    
    boxplot(data[[metric_col]] ~ data$method,
            col = "white",
            # Reuses the standard colors for the methods
            border = rep(2:7, times = 2), 
            ylab = y_label,
            xlab = "",              # Leave blank, names will act as x-axis labels
            las = 2,                # Rotates x-axis text to be perpendicular
            cex.axis = 0.85,        # Shrinks axis text slightly to fit
            ylim = c(y_min * 0.95, y_max * 1.05),
            main = "")
    
    # Add the dashed line for coverage plots (e.g., at 0.90)
    if (add_target_line) {
      abline(h = target_cov, lty = 2, col = "black", lwd = 1.5)
    }
  }

  # 3. Loop Through All Settings
  # ---------------------------------------------------------
  for (setting in setting_list) {
    cat(sprintf("Processing setting: %s...\n", setting))
    
    # Reset data frame for the current setting
    all_data <- data.frame()
    missing_files_count <- 0
    
    # Aggregate data from the folders
    for (i in 1:num_runs) {
      file_name <- sprintf("%s_seed_%d.csv", setting, i)
      file_path <- file.path(results_dir, as.character(i), file_name)
      
      # Check if the file exists before attempting to read it
      if (file.exists(file_path)) {
        temp_df <- read.csv(file_path, check.names = FALSE)
        all_data <- rbind(all_data, temp_df)
      } else {
        # Output a clean warning to the terminal if the file is missing
        cat(sprintf("  -> [WARNING] Missing file skipped: %s\n", file_path))
        missing_files_count <- missing_files_count + 1
      }
    }
    
    # Check if we actually found data for this setting
    if (nrow(all_data) == 0) {
      cat(sprintf("  -> [ERROR] No data found at all for %s. Skipping plot generation.\n", setting))
      next
    }
    
    if (missing_files_count > 0) {
      cat(sprintf("  -> Plotting %s using the %d available files...\n", setting, num_runs - missing_files_count))
    }
    
    # Ensure the methods plot in the correct order for the new data format
    method_order <- c("New CAMS", "CAMS", "DFT-adaptive-T (Joint)", "DFT-adaptive-CT (Joint)", "DFT-fixed (Joint)", 
                      "Vanilla CQR (Joint)", "Cox (Joint)",
                      "DFT-adaptive-T (Subgroup)", "DFT-adaptive-CT (Subgroup)", "DFT-fixed (Subgroup)", 
                      "Vanilla CQR (Subgroup)", "Cox (Subgroup)")
    all_data$method <- factor(all_data$method, levels = method_order)
    
    # 4. Generate the 6 Plots for the Current Setting as a PNG Image
    # ---------------------------------------------------------
    png_name <- file.path(plots_dir, sprintf("%s_results_plots.png", setting))
    
    # Open high-resolution PNG device (300 DPI is standard for publication)
    png(png_name, width = 15, height = 11, units = "in", res = 300) 
    
    # Set up a 2x3 grid. 
    # Added 'oma = c(0, 0, 3, 0)' to create an Outer Margin Area at the top for the global header
    par(mfrow = c(2, 3), mar = c(11, 4, 2, 1), oma = c(0, 0, 3, 0))
    
    # Plot 1: Marginal Coverage
    draw_paper_boxplot(all_data, 
                      metric_col = "Marginal coverage", 
                      y_label = "Overall Coverage Rate", 
                      add_target_line = TRUE)
    
    # Plot 2: Subgroup Coverage (X1 = 0)
    draw_paper_boxplot(all_data, 
                      metric_col = "group coverage for x_1 = 0", 
                      y_label = "Coverage rate (x_1 = 0)", 
                      add_target_line = TRUE)
    
    # Plot 3: Subgroup Coverage (X1 = 1)
    draw_paper_boxplot(all_data, 
                      metric_col = "group coverage for x_1 = 1", 
                      y_label = "Coverage rate (x_1 = 1)", 
                      add_target_line = TRUE)
    
    # Plot 4: Overall Average Lower Bound
    draw_paper_boxplot(all_data, 
                      metric_col = "lower bound values mean", 
                      y_label = "Overall Average Lower Bound", 
                      add_target_line = FALSE)

    # Plot 5: Lower Bound Mean for X1 = 0
    draw_paper_boxplot(all_data, 
                      metric_col = "lower bound mean for x_1 = 0", 
                      y_label = "Average Lower Bound (x_1 = 0)", 
                      add_target_line = FALSE)
                      
    # Plot 6: Lower Bound Mean for X1 = 1
    draw_paper_boxplot(all_data, 
                      metric_col = "lower bound mean for x_1 = 1", 
                      y_label = "Average Lower Bound (x_1 = 1)", 
                      add_target_line = FALSE)

    # Add the Global Header to the outer margin
    # 'outer = TRUE' places it in the 'oma' space, 'cex = 1.8' scales the text up, 'font = 2' makes it bold
    mtext(sprintf("Results for Setting: %s", setting), side = 3, outer = TRUE, cex = 1.8, font = 2)

    dev.off() # Close the PNG writer
    cat(sprintf("  -> Plot saved to: %s\n", png_name))
  }

  cat("All settings processed successfully!\n")
}

make_plots()