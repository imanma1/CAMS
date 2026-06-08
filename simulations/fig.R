# 1. Define Parameters
# ---------------------------------------------------------
results_dir <- "../results"
plots_dir <- file.path(results_dir, "plots")
num_runs <- 5             # The number of seed folders to loop through
target_alpha <- 0.1          # Miscoverage target
target_cov <- 1 - target_alpha

# List of all settings used in the simulations
setting_list <- c("ld_setting1", "ld_setting2", "ld_setting3", "ld_setting4",
                  "hd_homosc", "hd_heterosc")

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
          # Reuses the 6 standard colors for both Joint and Subgroup versions
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
  
  # Aggregate data from the folders
  for (i in 1:num_runs) {
    file_name <- sprintf("%s_seed_%d.csv", setting, i)
    file_path <- file.path(results_dir, as.character(i), file_name)
    
    if (file.exists(file_path)) {
      temp_df <- read.csv(file_path, check.names = FALSE)
      all_data <- rbind(all_data, temp_df)
    }
  }
  
  # Check if we actually found data for this setting
  if (nrow(all_data) == 0) {
    cat(sprintf("  -> No data found for %s, skipping.\n", setting))
    next
  }
  
  # Ensure the methods plot in the correct order for the new data format
  method_order <- c("CAMS", "DFT-adaptive-T (Joint)", "DFT-adaptive-CT (Joint)", "DFT-fixed (Joint)", 
                    "Vanilla CQR (Joint)", "Cox (Joint)", "Random Forest (Joint)",
                    "DFT-adaptive-T (Subgroup)", "DFT-adaptive-CT (Subgroup)", "DFT-fixed (Subgroup)", 
                    "Vanilla CQR (Subgroup)", "Cox (Subgroup)", "Random Forest (Subgroup)")
  all_data$method <- factor(all_data$method, levels = method_order)
  
  # 4. Generate the 6 Plots for the Current Setting
  # ---------------------------------------------------------
  pdf_name <- file.path(plots_dir, sprintf("%s_results_plots.pdf", setting))
  pdf(pdf_name, width = 15, height = 11) # Made PDF wider to comfortably fit 3 columns
  
  # Set up a 2x3 grid for the 6 plots, expanding bottom margin for longer text
  par(mfrow = c(2, 3), mar = c(11, 4, 2, 1))
  
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

  dev.off() # Close the PDF writer
  cat(sprintf("  -> Plots saved to: %s\n", pdf_name))
}

cat("All settings processed successfully!\n")