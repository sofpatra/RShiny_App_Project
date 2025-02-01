library('tidyverse')
library('biomaRt')
library('dplyr')
library('fgsea')


#create GMT file from drosophila GO file
# Read the GO annotation file
go_data <- read.delim("Data/Thermal/gene_association.fb", header = FALSE, comment.char = "!")



# Extract necessary columns: Gene ID (FBgn), GO Term
go_data <- go_data[, c(2, 5)]  # Column 2 = FBgn, Column 5 = GO Term
colnames(go_data) <- c("Gene", "GO")

# Group by GO term to create gene sets
gmt <- go_data %>%
  group_by(GO) %>%
  summarise(Genes = paste(Gene, collapse = "\t"))

# Format into GMT structure
gmt$Description <- "GO Annotation"  # Add a description 
gmt <- gmt[, c("GO", "Description", "Genes")]

# Save as a GMT file
write.table(gmt, "Data/Thermal/go_annotations.gmt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#Load in Deseq results

deseq_results <- read.csv('Data/Thermal/fluctuating_vs_control.csv')  # Read the CSV file


#Rank by Log Fold Change 
make_ranked_log2fc <- function(deseq_results) {
  # Create a tibble with log2 fold change and names
  log2fc_tibble <- tibble(GeneID = deseq_results$gene_id, log2FoldChange = deseq_results$log2FoldChange)
  
  # Get unique rows based on both GeneID and log2FoldChange
  unique_log2fc_tibble <- log2fc_tibble %>%
    distinct(GeneID, log2FoldChange) %>%  
    group_by(GeneID) %>%
    filter(n() == 1) %>%  # Keep only symbols that have a single unique log2FoldChange
    ungroup() %>%
    distinct(log2FoldChange, .keep_all = TRUE)
  
  # Convert to a named vector
  named_unique_log2fc <- setNames(unique_log2fc_tibble$log2FoldChange, unique_log2fc_tibble$GeneID)
  
  return(named_unique_log2fc)
}

rank <- make_ranked_log2fc(deseq_results)


#run GSEA
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- gmtPathways(gmt_file_path)
  
  fgseaRes <- tibble(fgsea(pathways, rnk_list, min_size, max_size))
  
  
  return(fgseaRes)
}

fgseaRes <- run_fgsea("Data/Thermal/go_annotations.gmt", rank, 10, 500)

# Connect to the Ensembl Mart
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl")

# Query GO term descriptions
go_descriptions <- getBM(attributes = c("go_id", "name_1006"), 
                         filters = "go", 
                         values = unique(fgseaRes$pathway), 
                         mart = ensembl)

# View the first few descriptions
head(go_descriptions)

simplified_fgsea_results <- fgseaRes %>%
  dplyr::select(-leadingEdge)  # Remove the 'leadingEdge' column

write.csv(simplified_fgsea_results, "Data/Thermal/gsea_results.csv", row.names = FALSE)
