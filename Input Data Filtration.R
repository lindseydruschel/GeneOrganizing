# Clear the environment
rm(list = ls())

# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(purrr) # For the reduce function

# Set file path (update with your actual file path)
file_path <- "C:/Users/druschel/Downloads/Mouse_GBM_full.xlsx"

# Get sheet names from the Excel file
sheet_names <- excel_sheets(file_path)

# Initialize an empty list to store data from sheets
data_list <- list()

# Loop through sheets 2 to the end
for (sheet in sheet_names[-1]) {
  # Read the current sheet
  sheet_data <- read_excel(file_path, sheet = sheet)
  
  # Check if the sheet contains valid data
  if (nrow(sheet_data) > 0) {
    # Select the desired columns and rename them to include the sheet name
    extracted_data <- sheet_data %>%
      select(gene, ave.norm.expr.1) %>%
      mutate(gene = toupper(gene)) %>% # Convert gene names to uppercase
      rename(
        !!paste0(sheet, "_ave.norm.expr.1") := ave.norm.expr.1
      )
    
    # Append the extracted data to the list
    data_list[[sheet]] <- extracted_data
  } else {
    message(paste("Sheet", sheet, "is empty. Skipping."))
  }
}

# Check if data_list is empty
if (length(data_list) == 0) {
  stop("No valid data found in any sheets.")
}

# Merge all data frames by the "gene" column
final_data <- purrr::reduce(data_list, full_join, by = "gene")

# Replace NAs with zeros
final_data[is.na(final_data)] <- 0

# View the combined data
print(head(final_data))

# Optionally, write the final data to a new CSV file
write.csv(final_data, "C:/Users/druschel/Downloads/Final_Combined_Data.csv", row.names = FALSE)
# Save the final data to a tab-delimited text file
write.table(final_data, "C:/Users/druschel/Downloads/Final_Combined_Data.txt", sep = "\t", row.names = FALSE, quote = FALSE)