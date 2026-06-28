merge <- function() {
  dir_final_results <- "../results"
  dir_results <- "../old_results"
  dir_results0 <- "../new_results"

  # Iterate through the subfolders 1 to 50
  for (folder_num in 1:50) {
    
    # Construct paths for the current subfolder
    current_res_dir <- file.path(dir_results, as.character(folder_num))
    current_res0_dir <- file.path(dir_results0, as.character(folder_num))
    current_final_dir <- file.path(dir_final_results, as.character(folder_num))
    
    # List all CSV files in the current results folder
    csv_files <- list.files(current_res_dir, pattern = "\\.csv$", full.names = FALSE)
    
    for (file_name in csv_files) {
      # Construct full file paths
      path_res <- file.path(current_res_dir, file_name)
      path_res0 <- file.path(current_res0_dir, file_name)
      path_final_res <- file.path(current_final_dir, file_name)
      
      # Check if the corresponding file actually exists in results0
      if (file.exists(path_res0)) {
        
        # Read both CSV files
        df_res <- read.csv(path_res, check.names = FALSE, stringsAsFactors = FALSE)
        df_res0 <- read.csv(path_res0, check.names = FALSE, stringsAsFactors = FALSE)
        
        # Extract the source row from results0
        cams_row <- df_res0[df_res0$method == "CAMS", ]
        
        # Proceed only if we found a source CAMS row to copy from
        if (nrow(cams_row) > 0) {
          
          # Rename the incoming row to "New CAMS" before inserting it
          cams_row$method <- "New CAMS"
          
          # Check if "New CAMS" already exists in the target file
          if ("New CAMS" %in% df_res$method) {
            
            # 1. Locate and remove the original "CAMS" row entirely
            old_cams_idx <- which(df_res$method == "CAMS")
            if (length(old_cams_idx) > 0) {
              df_res <- df_res[-old_cams_idx, ]
            }
            
            # 2. Rename the existing "New CAMS" row to "CAMS"
            # We do this by name rather than index so it works perfectly even after the deletion above
            df_res$method[df_res$method == "New CAMS"] <- "CAMS"
          }
          
          # 3. Bind the fresh incoming row to the VERY TOP of the dataframe
          # (This executes regardless of whether the if-statement triggered)
          df_res <- rbind(cams_row, df_res)
          
          # Re-index the first column so the numbering remains sequential and clean
          df_res[[1]] <- 1:nrow(df_res)
          
          dir.create(current_final_dir, showWarnings = FALSE, recursive = TRUE)
          write.csv(df_res, path_final_res, row.names = FALSE)
        }
      }
    }
  }

  print("Merge complete.")
}

switch_rows <- function() {
  dir_final_results <- "../results"

  for (folder_num in 1:50) {
    current_final_dir <- file.path(dir_final_results, as.character(folder_num))
    
    # List all CSV files in the current results folder
    csv_files <- list.files(current_final_dir, pattern = "\\.csv$", full.names = FALSE)
    for (file_name in csv_files) {
      # Construct full file paths
      path_final_res <- file.path(current_final_dir, file_name)

      if (file.exists(path_final_res)) {
        df <- read.csv(path_final_res, check.names = FALSE, stringsAsFactors = FALSE)
        # rownames(df)[rownames(df) == "CAMS"] <- "CAMStmp"
        # rownames(df)[rownames(df) == "New CAMS"] <- "CAMS"
        # rownames(df)[rownames(df) == "CAMStmp"] <- "New CAMS"

        df$method[df$method == "CAMS"] <- "CAMStmp"
        df$method[df$method == "New CAMS"] <- "CAMS"
        df$method[df$method == "CAMStmp"] <- "New CAMS"

        # row_order <- 1:nrow(df)
        # row_order[c(1, 2)] <- c(2, 1)
        # df_swapped <- df[row_order, ]

        write.csv(df, path_final_res, row.names = FALSE)
      }
    }
  }
}

copy_file <- function() {
  dir_results1 <- "../results0.1"
  dir_results3 <- "../results0.3"
  dir_cams_results1 <- "../cams_results0.1"
  dir_cams_results3 <- "../cams_results0.3"
  for (folder_num in 1:50) {
    dir1 <- file.path(dir_results1, as.character(folder_num))
    dir3 <- file.path(dir_results3, as.character(folder_num))
    dir_cams1 <- file.path(dir_cams_results1, as.character(folder_num))
    dir_cams3 <- file.path(dir_cams_results3, as.character(folder_num))

    file1 <- file.path(dir_cams1, sprintf("starve_hetero_high_dim_seed_%d.csv", folder_num))
    file3 <- file.path(dir_cams3, sprintf("starve_hetero_high_dim_seed_%d.csv", folder_num))

    file.copy(from = file1, to = dir1)
    file.copy(from = file3, to = dir3)
  }
}