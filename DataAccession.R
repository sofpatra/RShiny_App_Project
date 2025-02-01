library(tidyverse)
library(data.table)
library(GEOquery)

thermaldata <- fread("Data/Thermal/GSE150450_gene_count_matrix.csv.gz")

write.csv(thermaldata, "Data/Thermal/counts_data.csv", row.names = FALSE)

geo_data <- getGEO("GSE150450", GSEMatrix = TRUE)

geo_data

sample_metadata <- pData(geo_data[[1]])

sample_metadata <- sample_metadata %>%
  select(geo_accession, title, `lifestage:ch1`, `Sex:ch1`, `treatment:ch1`)

colnames(sample_metadata) <- c("geo_accession", "Sample", "Lifestage", "Sex", "Condition")
write.csv(sample_metadata, "Data/Thermal/metadata.csv", row.names = FALSE)
