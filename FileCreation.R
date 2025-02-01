
#load libraries
library(tidyverse)
library(data.table)

#read in files 
#create list of file names
files <- list.files(path = "Data/Prenatal_stress", pattern = "\\.gz$", full.names = TRUE)
#check
print(files)

merged_data <- NULL

#read the files into the list 
for (file in files) {
  # Read in the file
  data <- fread(file) %>%
    select(,c(1,7)) 
  # Extract the custom name by removing the path and extracting the last part after '_'
  base_name <- basename(file)
  custom_name <- sub(".*_(.*_.*)\\.gene\\.txt\\.gz$", "\\1", base_name)
  # Assign the data frame with the custom name
  colnames(data)[2] <- custom_name
  
  if (is.null(merged_data)) {
    merged_data <- data  # Initialize with the first file's data
  } else {
    merged_data <- left_join(merged_data, data, by = "Geneid")
  }
}

write.csv(merged_data, "Data/counts_data.csv", row.names = FALSE)

sample_names <- colnames(merged_data[,-1])

metadata <- data.frame(
  Sample = sample_names,
  Condition = ifelse(grepl("^Co", sample_names), "Control", "Stressed")
)

#Extract sex
metadata$Sex <- ifelse(grepl("_F", sample_names), "Female", "Male")

#View the generated metadata data frame
print(metadata)

write.csv(metadata, "Data/metadata.csv", row.names = FALSE)
