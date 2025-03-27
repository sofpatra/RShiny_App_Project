## Author: Sofiya Patra
## spatra@bu.edu
## BU BF591
## Project 

#Load in libraries
library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(gplots)
library(stats)
library(colourpicker) 
library(DT)


options(shiny.maxRequestSize=30*1024^2)

sample_columns <- c("Condition", "Sex", "Lifestage")
# Define UI for application that draws a histogram

columns = c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')
ui <- fluidPage(
  titlePanel("RNA-seq Analysis Application - Sofiya Patra"),
  tabsetPanel(
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(fileInput("metaupload", "Upload a metadata file", accept = c(".csv"))),
               mainPanel( tabsetPanel(
                          tabPanel("Summary",
                                  tableOutput("meta_summarytable")),
                          tabPanel("Table",
                                   DT::dataTableOutput("sampletable")),
                          tabPanel("Plots",
                                   radioButtons("x", "Choose sample variable to plot", sample_columns, selected = "Sex"),
                                   radioButtons("group", "Choose grouping variable", sample_columns, selected = "Condition"),
                                   submitButton(text= "Plot!", icon = icon("paint-brush"), width=400),
                                   plotOutput("sampleplot")))
             ))),
    tabPanel("Counts",
             sidebarLayout(
               sidebarPanel(fileInput("countsupload", "Upload a counts file", accept = c(".csv")),
                            "Filtering Parameters",
                            sliderInput("varianceslider", "Include genes with at least X percentile of variance", 0, 100, 20),
                            sliderInput("zeroslider","Include genes with at least X samples that are non-zero", 0, 80, 40),
                            submitButton(text= "Filter!", icon = icon("paint-brush"), width=550)),
               mainPanel( tabsetPanel(
                          tabPanel("Filtering",
                                   sidebarLayout(
                                     sidebarPanel("Summary Table",
                                                  tableOutput("counts_summarytable")),
                                     mainPanel("Counts Table",
                                               tableOutput("countstable")))),
                          tabPanel("Diagnostic Plots",
                                   plotOutput("scatter1"),
                                   plotOutput("scatter2")),
                          tabPanel("Heatmap",
                                   plotOutput("heatmap")),
                          tabPanel("PCA",
                                   radioButtons("shape", "Choose shape mapping", sample_columns, selected = "Sex"),
                                   radioButtons("pcacolor", "Choose color mapping", sample_columns, selected = "Lifestage"),
                                   submitButton(text= "Plot!", icon = icon("paint-brush"), width=400),
                                   plotOutput("pcaplot")))))),
    tabPanel("Differential Expression",
             sidebarLayout(
               sidebarPanel(fileInput("desequpload", "Upload a differential expression results file", accept = c(".csv")),
                            'A volcano plot can be generated with "log2 fold-change" on the x-axis and "p-adjusted" on the y-axis.',
                            radioButtons("xaxis", "Choose the column for the x axis", columns, selected = "log2FoldChange"),
                            radioButtons("yaxis", "Choose the column for the y axis", columns, selected = "padj"),
                            colourInput("color1", "Base point color", value = '#FFE203'),
                            colourInput("color2", "Highlight point color", value = '#347BD9'), 
                            sliderInput("slider", "Select the magnitude of the p adjusted coloring:", -30, 0, -10),
                            submitButton(text= "Plot!", icon = icon("paint-brush"), width=400)),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Plot",
                            plotOutput("volcano")),
                   tabPanel("Table",
                            DT::dataTableOutput("deseqtable"))
                 ))
             )),
    tabPanel("Gene Set Enrichment Analysis",
             tabsetPanel(
               tabPanel("Barplot",
                        sidebarLayout(
                          sidebarPanel(fileInput("gsea", "Upload a GSEA results table", accept =c(".csv")),
                                       sliderInput("top_pathways", "Number of pathways to include:", 1, 20, 5),
                                       submitButton(text= "Plot!", icon = icon("paint-brush"), width=400)),
                          mainPanel("Plot",
                                    plotOutput("gsea_barplot")))),
               tabPanel("GSEA Table",
                        sidebarLayout(
                          sidebarPanel(sliderInput("gseaslider2", "Select adjusted p-value filter:", -30, 0, -3),
                                       radioButtons("gsea_buttons", "Select pathways:", c("positive", "negative", "all")),
                                       submitButton(text= "Show Me!", icon = icon("paint-brush"), width=400),
                                       downloadButton("download_gsea", "Download Filtered Data")),
                          mainPanel(DT::dataTableOutput("gsea_table1"))
                        )),
               tabPanel("Scatter Plot",
                        sidebarLayout(
                          sidebarPanel(sliderInput("gseaslider3", "Select adjusted p-value filter:", -30, 0, -5),
                                       submitButton(text= "Plot!", icon = icon("paint-brush"), width=400)),
                          mainPanel(plotOutput("gsea_scatterplot"))
                        )))
             )
  )
  
  
)


# Define server logic 
server <- function(input, output, session) {
  
  #load in metadata file
  load_metadata <- reactive({
    req(input$metaupload)  # Ensure file input is available
    metadata <- read.csv(input$metaupload$datapath, stringsAsFactors = TRUE)  # Read the CSV file
  })
  
  # load in counts data file
  load_countsdata <- reactive({
    req(input$countsupload)  # Ensure file input is available
    countsdata <- read.csv(input$countsupload$datapath)  # Read the CSV file
  })
  
  #load in deseq results file
  load_deseq <- reactive({
    req(input$desequpload)  # Ensure file input is available
    countsdata <- read.csv(input$desequpload$datapath)  # Read the CSV file
  })
  
  #load in gsea results table 
  load_gsea <- reactive({
    req(input$gsea)
    gsea_data <- read.csv(input$gsea$datapath)
  })
  
  #create a sample summary table 
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
      })
    )
    return(summary)
  }
  
  sample_plot <-
    function(metadata, x_name, group_name) {
      plot <-  ggplot(metadata, aes_string(x= x_name, fill = group_name)) +
        geom_bar(position = "dodge") +
        labs(x= x_name,y= "Count", title = "Sample Data")+
        theme_minimal()
      return(plot)
    }
  
  # Filtering 
  filtered_countsdata <- reactive({
    countsdata <- load_countsdata()
    req(countsdata)
    #calculate median for each gene 
    countsdata$Median <- apply(countsdata[, -1], 1, median)
    #Calculate variance for each gene
    countsdata$Variance <- apply(countsdata[, -1], 1, var)
    #Count Zeros
    countsdata$Zeros <- rowSums(countsdata[, -1] == 0)
    #Filter
    filtered <- countsdata[
      countsdata$Variance >= quantile(countsdata$Variance, probs = input$varianceslider /100) &
        (ncol(countsdata) - countsdata$Zeros -1) >= input$zeroslider, ]
    
    return(list(filtered = filtered, all_data = countsdata))
  })
  
  filter_summary <- reactive({
    data <- filtered_countsdata()
    filtered <- data$filtered
    all_data <- data$all_data
    
    # Calculate summary statistics
    num_samples <- ncol(all_data[2:81])
    num_genes_total <- nrow(all_data)
    num_genes_filtered <- nrow(filtered)
    percent_filtered <- round(100 * num_genes_filtered / num_genes_total, 2)
    
    summary <- data.frame(
      Metric = c("Number of Samples", "Total Genes", "Genes Passing Filters", "Genes Failing Filters"),
      Count = c(num_samples, num_genes_total, num_genes_filtered, num_genes_total - num_genes_filtered),
      Percentage = c(100, 100, percent_filtered, 100 - percent_filtered)
    )
    
    return(summary)
  })
  
  #Filtering Plots 
  variance_plot <-
    function(data) {
      plot <-  ggplot(data, aes(x= Median, y = Variance, color = Pass)) +
        geom_point(alpha = 0.8) +
        scale_color_manual(values = c("Pass" = "darkblue", "Fail" = "gray")) +
        labs(x = "Median Count (Log Scale)", y = "Variance", title = "Median Count vs Variance") +
        theme_minimal() + 
        scale_x_log10() +
        scale_y_log10()
      
      return(plot)
    }
  
  zeros_plot <- 
    function(data) {
      plot <- ggplot(data, aes(x= Median, y = Zeros, color = Pass)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Pass" = "darkblue", "Fail" = "gray")) +
      labs(x = "Median Count (Log Scale)", y = "Number of Zeros", title = "Median Count vs Number of Zeros") +
      theme_minimal() +
      scale_x_log10()
      
      return(plot)
    }
  
  heatmap <- 
    function(data) {
      mat <- as.matrix(data[, 2:81])
      rownames(mat) <- rep("", nrow(mat))
      heatmap.2(mat, scale = "row",
                col="bluered",  # colors
                labRow = NULL,
                srtCol = 90,            # Rotation angle for column labels (0 for horizontal)
                cexCol = 0.7,         # Column label size
                trace = "none",       # Removes trace lines in the heatmap
                dendrogram = "column"
      )
    }
  
  plot_pca <- 
    function(data, metadata, shape_column, color_column) {
      expr_mat <- data[, !colnames(data) %in% c("Median", "Zeros", "Variance")] %>%
        pivot_longer(-c(gene_id), names_to = "sample") %>%
        pivot_wider(names_from = gene_id)
      
      
      expr_mat <- expr_mat %>%
        column_to_rownames(var = "sample") %>%
        as.matrix()
      #perform PC
      print("Starting PCA...")
      pca <- prcomp(
        expr_mat,
        center=TRUE,
        scale=TRUE
      )
      print("PCA completed successfully!")
      
      # Calculate variance explained
      variance_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
      
      # Format axis labels with variance explained
      pc1_label <- paste0("PC1 (", round(variance_explained[1], 2), "%)")
      pc2_label <- paste0("PC2 (", round(variance_explained[2], 2), "%)")
      
      
      pca_results <- as_tibble(pca$x, rownames = "Sample") %>%
        left_join(metadata, by = "Sample")  # Join metadata
      
      plot <- pca_results %>%
        ggplot(aes(x=PC1,y=PC2, color = .data[[color_column]], shape = .data[[shape_column]])) +
        geom_point() +
        labs(title = "PCA Plot", x= pc1_label, y= pc2_label) +
        theme_minimal()
      
      return(plot)
    }
  #filter GSEA
  filtered_gsea <- reactive({
    req(load_gsea())
    gsea_data <- load_gsea()
    
    gsea_data <- gsea_data %>% filter(!is.na(padj))
    #order by padj 
    gsea_data <- gsea_data %>%
      arrange(padj) %>%
      head(input$top_pathways)
    
    return(gsea_data)
  })
  
  #Second Tab Filter
  filtered_gsea2 <- reactive({
    req(load_gsea())
    gsea_data <- load_gsea()
    #P-value filter
    data_filtered <- gsea_data %>%
      filter(padj <= 10^input$gseaslider2) %>%
      mutate(
        pval = formatC(pval, format ="e", digits = 5),
        padj = formatC(padj, format = "e", digits = 5)
      ) %>%
      arrange(padj)
    # Apply NES filter
    if (input$gsea_buttons == "positive") {
      data_filtered <- data_filtered %>% filter(NES > 0)
    } else if (input$gsea_buttons == "negative") {
      data_filtered <- data_filtered %>% filter(NES < 0)
    }
    return(data_filtered)
  })
  
  #create gsea barplot 
  bar_plot <- 
    function(data) {
      ggplot(data, aes(x = reorder(pathway, NES), y= NES, fill = NES)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = "Top Pathways by NES",
             x= "Pathway",
             y= "Normalized Enrichment Score (NES)") +
        theme_minimal()
    }
  
  #' Volcano plot
  #'
  #' @param dataf The loaded data frame.
  #' @param x_name The column name to plot on the x-axis
  #' @param y_name The column name to plot on the y-axis
  #' @param slider A negative integer value representing the magnitude of
  #' p-adjusted values to color. Most of our data will be between -1 and -300.
  #' @param color1 One of the colors for the points.
  #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
      plot <-  ggplot(dataf, aes_string(x= x_name, y=paste0("-log10(", y_name, ")"))) +
        geom_point(aes(color = !!sym(y_name) < 10^slider)) +
        scale_color_manual(values = c("TRUE" = color2, "FALSE" = color1)) +
        labs(x= x_name, y =  paste("-log10(", y_name, ")"), title = "Volcano Plot")+
        theme_minimal()
      
      
      return(plot)
    }
  
  #' Draw and filter table
  #' @param dataf Data frame loaded by load_data()
  #' @param slider Negative number, typically from the slider input.
  #'
  #' @return Data frame filtered to p-adjusted values that are less than 
  #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
  #' displayed.
  draw_table <- function(dataf, slider) {
    filtered <- dataf %>%
      filter(padj < 1 * 10^(slider)) %>%
      mutate(
        pvalue = formatC(pvalue, format ="e", digits = 5),
        padj = formatC(padj, format = "e", digits = 5)
      )
    return(filtered)
  }
  
  #Outputs 

  # Output for sample table 
  output$sampletable <- DT::renderDataTable({
    metadata <- load_metadata()
    req(metadata)
    DT::datatable(metadata, options = list(pageLength = 10))
  })
  
  
  # Same here, just return the table as you want to see it in the web page
  output$meta_summarytable <- renderTable({
    metadata <- load_metadata()
    req(metadata)
    create_summary_table(metadata)

  }) 
  
  output$sampleplot <- renderPlot({
    metadata <- load_metadata()
    req(metadata)
    sample_plot(metadata, input$x, input$group)
  })
  
  #Outputs for counts  
  output$countstable <- renderTable({
    filteredcounts <- filtered_countsdata()$filtered
    req(filteredcounts)
    head(filteredcounts, 50)
  })
  
  output$counts_summarytable <- renderTable({
    req(filter_summary()) 
    filter_summary() 
  })
  
  output$scatter1 <- renderPlot({
    data <- filtered_countsdata()$all_data
    data$Pass <- ifelse(data$Variance >= quantile(data$Variance, probs = input$varianceslider / 100) &
                          (ncol(data) - data$Zeros -1) >= input$zeroslider, "Pass", "Fail")
    variance_plot(data)
  })
  
  output$scatter2 <- renderPlot({
    data <- filtered_countsdata()$all_data
    data$Pass <- ifelse(data$Variance >= quantile(data$Variance, probs = input$varianceslider / 100) &
                          (ncol(data) - data$Zeros -1) >= input$zeroslider, "Pass", "Fail")
    zeros_plot(data)
  }) 
  
  output$heatmap <- renderPlot({
    data <- filtered_countsdata()$filtered
    heatmap(data)
  })
  
  output$pcaplot <- renderPlot({
    data <- filtered_countsdata()$filtered
    metadata <- load_metadata()
    req(metadata)
    req(data)
    shape_column <- input$shape
    color_column <- input$pcacolor
    plot_pca(data, metadata, shape_column, color_column)
  })
  
  #Deseqtable output
  output$deseqtable <- DT::renderDataTable({
    data <- load_deseq()
    req(data)
    DT::datatable(data, options = list(pageLength = 10))
  })
  #Deseq Volcano Plot output
  output$volcano <- renderPlot({
    data <- load_deseq()
    req(data)
    volcano_plot(data, input$xaxis, input$yaxis, input$slider, color1 = input$color1, color2 = input$color2)
  }) 
  
  #GSEA Barplot output
  output$gsea_barplot <- renderPlot({
    req(filtered_gsea())
    gsea_data <- filtered_gsea()
    bar_plot(gsea_data)
  })
  #GSEA Table output1
  output$gsea_table1 <- DT::renderDataTable({
    req(filtered_gsea2())
    datatable(filtered_gsea2(), options = list(pageLength = 10, order = list(list(1, 'asc'))))
  })
  
  output$download_gsea <- downloadHandler(
    filename = function() {
      paste("filtered_gsea_results.csv")
    },
    content = function(file) {
      write.csv(filtered_gsea2(), file, row.names = FALSE)
    }
  )
  
  output$gsea_scatterplot <- renderPlot({
    req(load_gsea())
    data <- load_gsea()
    
    # Apply the p-value filter
    data <- data %>%
      filter(padj <= 10^input$gseaslider3)
    
    # Plot the scatter plot of NES vs. -log10 p-value
    ggplot(data, aes(x = NES, y = -log10(padj))) +
      geom_point(aes(color = NES), size = 3) +
      labs(title = "NES vs. -log10 Adjusted P-value", x = "NES", y = "-log10(Adjusted P-value)") +
      theme_minimal()
  })
  
}


# Run the application
shinyApp(ui = ui, server = server)