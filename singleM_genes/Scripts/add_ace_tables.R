# Set the folder path
folder_path <- "singleM_genes/Data/diversity_out/"

# Get a list of files with names containing ".chao1.csv"
files <- list.files(folder_path, pattern = ".chao1.csv", full.names = TRUE)

# Loop through the files
for (file in files) {
  # Read the original table
  original_table <- read.csv(file)
  
  # Duplicate the table (assuming it's a data frame)
  duplicated_table <- original_table
  
  # Rename columns or make any modifications as needed
  # For example, renaming ".chao1.csv" to ".ace.csv"
  new_filename <- gsub(".chao1.csv", ".ace.csv", basename(file))
  
  # Save the duplicated and renamed table to the same folder
  write.csv(duplicated_table, file.path(folder_path, new_filename), row.names = FALSE)
  
  # Print a message for each file processed
  cat("Duplicated and renamed:", basename(file), "->", new_filename, "\n")
}
