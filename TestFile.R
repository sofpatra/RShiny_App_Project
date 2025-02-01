library(bslib)
library(ggplot2)
library(tidyverse)


#Loading Data
metadata <- read.csv('Data/Thermal/metadata.csv', stringsAsFactors = TRUE)  # Read the CSV file


countsdata <- read.csv('Data/Thermal/counts_data.csv')  # Read the CSV file

countsdata2 <- read.csv('Data/counts_data2.csv')
#Creating Summary Table
create_summary_table <- function(data) {
  summary <- data.frame(
    Column_Name = names(data),
    Type = sapply(data, function(x) class(x)[1]),
    Summary = sapply(data, function(x) {
      if (is.numeric(x)) {
        paste0(round(mean(x, na.rm = TRUE), 2), " (", round(sd(x, na.rm = TRUE), 2), ")")
      } else {
        paste(unique(x), collapse = ", ")
      }
    }),
    stringsAsFactors = FALSE
  )
  return(summary)
}

create_summary_table(metadata)

#Create Sample Plot 
sample_plot <-
  function(metadata, x_name) {
    plot <-  ggplot(metadata, aes_string(x= x_name)) +
      geom_histogram(stat = "count") +
      labs(x= x_name, title = "Sample Data")+
      theme_minimal()
    return(plot)
  }

sample_plot(metadata, "Condition")


countsdata$Median <- apply(countsdata[, -1], 1, median)
#Calculate variance for each gene
countsdata$Variance <- apply(countsdata[, -1], 1, var)
#Count Zeros
countsdata$Zeros <- rowSums(countsdata[, -1] == 0)
#Filter
filtered <- countsdata[
    countsdata$Variance >= quantile(countsdata$Variance, probs = 20 /100) &
      (ncol(countsdata) - countsdata$Zeros -1) >= 6, ]
colnames(filtered)

#heatmap 
heatmap <- function(data) {
  mat <- as.matrix(data[, 2:81])
  rownames(mat) <- NULL
  print(rownames(mat))
  heatmap.2(mat, scale = "row",
            col="bluered",  # colors
            labRow = NULL,
            srtCol = 90,            # Rotation angle for column labels (0 for horizontal)
            cexCol = 0.7,         # Column label size
            trace = "none",       # Removes trace lines in the heatmap
            dendrogram = "column"
            )
}
heatmap(filtered)


#PCA
expr_mat <- filtered %>%
      select(-c(Median, Zeros,Variance)) %>%
      pivot_longer(-c(gene_id), names_to = "sample") %>%
      pivot_wider(names_from = gene_id)


expr_mat <- expr_mat %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()
    
pca <- prcomp(
      expr_mat,
      center=TRUE,
      scale=TRUE
    )

pca_results <- as_tibble(pca$x, rownames = "Sample") %>%
  left_join(metadata, by = "Sample")  # Join metadata

colnames(pca_results)

plot <- pca_results %>%
      ggplot(aes(x=PC1,y=PC2, color = Lifestage, shape = Sex)) +
      geom_point() +
      labs(title = "PCA Plot") +
      theme_minimal()
    
plot


#Deseq 
results <- read.csv('Data/Thermal/fluctuating_vs_control.csv')  # Read the CSV file
volcano_plot <-
  function(dataf, x_name, y_name, slider, color1, color2) {
    plot <-  ggplot(dataf, aes_string(x= x_name, y=paste0("-log10(", y_name, ")"))) +
      geom_point(aes(color = !!sym(y_name) < 10^slider)) +
      scale_color_manual(values = c("TRUE" = color2, "FALSE" = color1)) +
      labs(x= x_name, y =  paste("-log10(", y_name, ")"), title = "Volcano Plot")+
      theme_minimal()
    
    
    return(plot)
  }
volcano_plot(results, 'log2FoldChange', 'padj', -10, color1= 'blue', color2 = 'pink')


#GSEA 
gsea_results <- read.csv('Data/Thermal/gsea_results.csv')


