#load libraries
library(ggplot2)
library(tidyverse)
library(data.table)
library(DESeq2)
library(biomaRt)

metadata <- read.csv('Data/Thermal/metadata.csv', stringsAsFactors = TRUE)  # Read the CSV file


countsdata <- read.csv('Data/Thermal/counts_data.csv')  # Read the CSV file

filtered_counts <- countsdata[rowSums(countsdata[-1] != 0) >= 50, ]

count_mat <- as.matrix(filtered_counts[-1])
row.names(count_mat) <- filtered_counts$gene_id

metadata$Sex_Lifestage <- interaction(metadata$Sex, metadata$Lifestage, drop = TRUE)

design <- formula(~ Sex_Lifestage + Condition)

dds <- DESeqDataSetFromMatrix(
  countData=count_mat,
  colData=metadata,
  design=design
)
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name="Condition_fluctuating_vs_control")
fluctuating_vs_control_de <- as_tibble(res) %>%
  mutate(gene_id=rownames(res)) %>%
  relocate(gene_id) %>%
  arrange(pvalue)

# Connect to the Ensembl Mart
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl")

mart <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

flybase_ids <- res$gene_id 

id_mapping <- getBM(
  attributes = c("flybase_gene_id", "external_gene_name"),  # Columns to retrieve
  filters = "flybase_gene_id",                             # Filter type
  values = flybase_ids,                                    # Your list of FlyBase Gene IDs
  mart = mart
)

colnames(id_mapping) <- c("gene_id", "Symbol")


#add id column to results 
mapped_results <- res %>%
  left_join(id_mapping)


write.csv(mapped_results, "Data/Thermal/fluctuating_vs_control.csv", row.names = FALSE)
