# Load necessary libraries
if (!require(shiny)) install.packages("shiny")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(DT)) install.packages("DT")

source("preprocess_data.R")

# === Load Libraries ===
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(shinyWidgets)
library(clusterProfiler)

merged_sing_df <- readRDS("merged_sing_df.rds")
gmt_data <- clusterProfiler::read.gmt("data/20241021_188genelist_withphenotypes.gmt")
gmt_data
# === Define UI ===
ui <- fluidPage(
  titlePanel("Pathway Enrichment Comparison"),
  sidebarLayout(
    sidebarPanel(
      pickerInput(
        "pathway",
        "Select Pathways:",
        choices = unique(merged_sing_df$Pathway),
        selected = unique(merged_sing_df$Pathway)[1],
        options = list(`actions-box` = TRUE),
        multiple = TRUE
      ),
      selectInput("study", "Select Study", 
                  choices = c("All", unique(merged_sing_df$study))),
      selectInput("comparison", "Select Comparison Type", 
                  choices = c("Response Status" = "Response", "Recurrence Status" = "Recurrence Status")),
      selectInput("timepoint", "Select Timepoint", choices = c("Both", "Baseline", "Week 6")),
      downloadButton("downloadPlot", "Download plots"),
      downloadButton("downloadTable", "Download Table"),
      h4("Genes in Selected Signature(s)"),
      uiOutput("geneList"),
      width = 3
    ),
    mainPanel(
      fluidRow(
        column(12, DT::dataTableOutput("dataTable"))
      ),
      fluidRow(
        column(12, uiOutput("dynamicPlots"))
      )
    )
  )
)

# === Define Server Logic ===
server <- function(input, output, session) {
  
  # Reactive expression to get associated genes from GMT
  output$geneList <- renderUI({
    selected_pathways <- input$pathway
    
    if (length(selected_pathways) == 0) {
      return(div("No pathways selected."))
    }
    
    # Get genes for each selected pathway
    gene_lists <- lapply(selected_pathways, function(pathway) {
      genes <- gmt_data$gene[gmt_data$term == pathway]
      if (length(genes) > 0) {
        list(
          tags$h5(strong(pathway)), # Header for the pathway
          tags$ul(lapply(genes, tags$li)) # List of genes
        )
      } else {
        list(
          tags$h5(strong(pathway)),
          div("No genes found for this pathway.")
        )
      }
    })
    
    # Return combined UI
    do.call(tagList, gene_lists)
  })
  
  # Reactive data subset for the plot and table
  selected_data <- reactive({
    data <- merged_sing_df %>%
      filter(Pathway %in% input$pathway) %>%
      mutate(
        Comparison = case_when(
          input$comparison == "Response" ~ ifelse(MPRvNMPR == 1, "MPRs", "NMPRs"),
          input$comparison == "Recurrence Status" ~ ifelse(recurrence_status == 1, "Recurrence", "No Recurrence")
        )
      )
    
    # Subset by timepoint if not "Both"
    if (input$timepoint != "Both") {
      data <- data %>% filter(Timepoint == input$timepoint)
    }
    
    # Subset by study if not "All"
    if (input$study != "All") {
      data <- data %>% filter(study == input$study)
    }
    
    return(data)
  })
  
  # Generate Plot
  output$dynamicPlots <- renderUI({
    plot_list <- lapply(input$pathway, function(path) {
      plotname <- paste0("plot_", path)
      plotOutput(plotname, height = "600px", width = "100%")
    })
    do.call(tagList, plot_list)
  })
  
  # Render plots for each selected pathway
  observe({
    lapply(input$pathway, function(path) {
      output[[paste0("plot_", path)]] <- renderPlot({
        
        # Filter data for the selected pathway
        data <- selected_data() %>% filter(Pathway == path)
        
        # Compute Wilcoxon p-values for each study and timepoint
        p_values_df <- data %>%
          group_by(Timepoint, study) %>%
          summarise(
            p_value = tryCatch(
              wilcox.test(Singscore ~ Comparison)$p.value,
              error = function(e) NA
            ),
            .groups = 'drop'
          ) %>%
          mutate(
            label = case_when(
              is.na(p_value)          ~ "NA",
              p_value < 0.0001        ~ "p < 0.0001",
              p_value < 0.001         ~ "p < 0.001",
              p_value < 0.01          ~ "p < 0.01",
              p_value < 0.05          ~ "p < 0.05",
              TRUE                    ~ paste0("p = ", round(p_value, 3))
            )
          )
        
        # Plot
        ggplot(data, aes(x = Timepoint, y = Singscore, fill = Comparison)) +
          geom_boxplot(width = 0.5, alpha = 0.8, position = position_dodge(0.7), outlier.shape = NA) +
          geom_jitter(
            shape = 21, size = 3,
            position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
            color = "black", aes(fill = Comparison)
          ) +
          scale_fill_manual(values = if (input$comparison == "Response") {
            c("MPRs" = "#0072B2", "NMPRs" = "red")
          } else {
            c("No Recurrence" = "#56B4E9", "Recurrence" = "#D55E00")
          }) +
          facet_wrap(~study) +
          geom_text(
            data = p_values_df,
            aes(x = Timepoint, y = 0.6, label = label),
            inherit.aes = FALSE,
            size = 5
          ) +
          theme_minimal() +
          labs(
            title = paste(gsub("_", " ", path), "Enrichment"),
            x = "Time",
            y = paste(path, "Singscore"),
            fill = input$comparison
          ) +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            panel.spacing = unit(0.5, "lines"),
            strip.background = element_rect(fill = "grey80", color = "black"),
            strip.text = element_text(face = "bold", size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12)
          )
      })
    })
  })
  
  # Render Data Table
  output$dataTable <- DT::renderDataTable({
    DT::datatable(selected_data())
  })
  
  # Download all plots as a ZIP
  output$downloadPlot <- downloadHandler(
    filename = function() {
      "Pathway_Enrichment_Plots.zip"
    },
    content = function(file) {
      temp_dir <- tempdir()
      file_paths <- c()
      
      lapply(input$pathway, function(path) {
        plot_path <- file.path(temp_dir, paste0(path, ".png"))
        
        ggsave(filename = plot_path, plot = output[[paste0("plot_", path)]](), width = 10, height = 5, dpi = 300)
        
        file_paths <<- c(file_paths, plot_path)
      })
      
      zip::zip(zipfile = file, files = file_paths)
    },
    contentType = "application/zip"
  )
  
  # Download Table as CSV
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("Pathway_Enrichment_Data_", paste(input$pathway, collapse = "_"), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(selected_data(), file, row.names = FALSE)
    }
  )
}

# === Run the Application ===
shinyApp(ui = ui, server = server)