# Define the base directories
dir_results <- "../results"
dir_results0 <- "../results0"

# Iterate through the subfolders 1 to 5
for (folder_num in 1:50) {
  
  # Construct paths for the current subfolder
  current_res_dir <- file.path(dir_results, as.character(folder_num))
  current_res0_dir <- file.path(dir_results0, as.character(folder_num))
  
  # List all CSV files in the current results folder
  csv_files <- list.files(current_res_dir, pattern = "\\.csv$", full.names = FALSE)
  
  for (file_name in csv_files) {
    # Construct full file paths
    path_res <- file.path(current_res_dir, file_name)
    path_res0 <- file.path(current_res0_dir, file_name)
    
    # Check if the corresponding file actually exists in results0
    if (file.exists(path_res0)) {
      
      # Read both CSV files
      df_res <- read.csv(path_res, check.names = FALSE, stringsAsFactors = FALSE)
      df_res0 <- read.csv(path_res0, check.names = FALSE, stringsAsFactors = FALSE)
      
      # Extract the replacement row from results0
      cams_row <- df_res0[df_res0$method == "CAMS", ]
      
      # Proceed only if we found a source CAMS row to copy from
      if (nrow(cams_row) > 0) {
        
        # Locate the row index of the existing "CAMS" row in the target file
        target_idx <- which(df_res$method == "CAMS")
        
        # Proceed only if the "CAMS" row actually exists in the target file to overwrite
        if (length(target_idx) > 0) {
          
          # Replace the row data in-place
          df_res[target_idx, ] <- cams_row
          
          # Re-index the first column so the numbering remains sequential and clean
          df_res[[1]] <- 1:nrow(df_res)
          
          # Overwrite the original file in the 'results' folder with the updated data
          write.csv(df_res, path_res, row.names = FALSE)
        }
      }
    }
  }
}

print("Replacement complete!")