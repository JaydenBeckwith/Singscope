# Load necessary libraries
if (!require(shiny)) install.packages("shiny")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(DT)) install.packages("DT")

options(shiny.maxRequestSize = 400*1024^2)  # 400 MB
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
  titlePanel("Signature Enrichment Comparison"),
  
  # Main Tab Layout
  tabsetPanel(
    tabPanel(
      "Data Import",
      fluidRow(
        column(
          6, offset = 3,  # Centers the panel horizontally
          wellPanel(
            h3("Import Data"),
            fileInput("exprMatrix", "Upload Gene Expression Matrix (.csv, .tsv)", 
                      accept = c(".csv", ".tsv")),
            fileInput("metadata", "Upload Metadata (.csv, .tsv)", 
                      accept = c(".csv", ".tsv")),
            textInput("cohortName", "Enter Cohort Name", value = "NewCohort"),
            actionButton("submitData", "Submit Data"),
            div(style = "margin-bottom: 20px;")

          )
        ),
       
      ),
      fluidRow(
        column(
          12,
          h3("Uploaded Data Preview"),
          DT::dataTableOutput("previewExprMatrix"),
          DT::dataTableOutput("previewMetadata")
        )
      )
    ),
    
    tabPanel(
      "Signature Analysis",
      sidebarLayout(
        sidebarPanel(
          # Other Sidebar Options
          pickerInput(
            "pathway",
            "Select Signatures",
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
            column(12, 
                   h3("Global Top Statistically Significant Comparisons"),
                   div(style = "margin-bottom: 20px;"),
                   DT::dataTableOutput("globalTopSignificantTable")
            )
          ),
          fluidRow(
            div(style = "height: 30px;")  # 30px of vertical space
          ),
          fluidRow(
            column(12,
                   h3("Sample Distribution"),
                   plotOutput("sampleDistributionPlot", height = "400px"),
                   div(style = "margin-bottom: 20px;")
            )
          ),
          fluidRow(
            div(style = "height: 30px;")  # 30px of vertical space
          ),
          fluidRow(
            column(12,
                   h3("Selected Sample Data"),
                   DT::dataTableOutput("dataTable")
            )
          ),
          fluidRow(
            h3("Selected Signature Plots"),
            column(12, uiOutput("dynamicPlots"))
          )
        )
      )
    )
  )
)




# === Define Server Logic ===
server <- function(input, output, session) {
  

  # === Upload and Preview Data ===
  expr_data <- reactiveVal(NULL)
  metadata_data <- reactiveVal(NULL)
  
  observeEvent(input$submitData, {
    req(input$exprMatrix, input$metadata)
    
    # === 1. Load the Data ===
    expr <- read.csv(input$exprMatrix$datapath, row.names = 1)
    meta <- read.csv(input$metadata$datapath)
    
    # === 2. Save to Reactive Values for Display ===
    expr_data(expr)
    metadata_data(meta)
    
    # === 3. Update the UI with Upload Status ===
    output$uploadStatus <- renderUI({
      div(
        h4("Data Successfully Imported"),
        p(paste("Cohort Name: ", input$cohortName))
      )
    })
    
    output$previewExprMatrix <- DT::renderDataTable({
      DT::datatable(expr_data(), options = list(pageLength = 5))
    })
    
    output$previewMetadata <- DT::renderDataTable({
      DT::datatable(metadata_data(), options = list(pageLength = 5))
    })
    
    # === 4. Preprocessing Status Update ===
    output$preprocessStatus <- renderUI({
      div(
        h4("Preprocessing in progress..."),
        p("Please wait while your data is being processed.")
      )
    })
    
    tryCatch({
      # Call the preprocessing function and append to RDS
      new_data <- preprocess_data(
        exprMatrixPath = input$exprMatrix$datapath,
        metadataPath = input$metadata$datapath,
        cohortName = input$cohortName
      )
      
      # 📝 Update the global reference in memory
      merged_sing_df <<- new_data
      
      # Refresh Pathway and Study options
      updatePickerInput(session, "pathway", choices = unique(merged_sing_df$Pathway))
      updateSelectInput(session, "study", choices = c("All", unique(merged_sing_df$study)))
      
      # Show success notification
      showNotification("Data imported and merged successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  
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
  
  observeEvent(input$submitData, {
    req(input$exprMatrix, input$metadata)
    
    tryCatch({
      # Call the preprocessing function and append to RDS
      new_data <- preprocess_data(
        exprMatrixPath = input$exprMatrix$datapath,
        metadataPath = input$metadata$datapath,
        cohortName = input$cohortName
      )
      
      
      # Refresh Pathway and Study options
      updatePickerInput(session, "pathway", choices = unique(new_data$Pathway))
      updateSelectInput(session, "study", choices = c("All", unique(new_data$study)))
      
      # Show success notification
      showNotification("Data imported and merged successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
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
  
  global_comparisons <- reactive({
    # Perform the Wilcoxon test for both comparisons
    response_status <- merged_sing_df %>%
      group_by(Pathway, Timepoint, study) %>%
      summarise(
        `Comparison Type` = "Response Status",
        p_value = tryCatch(
          wilcox.test(Singscore ~ MPRvNMPR)$p.value,
          error = function(e) NA
        ),
        .groups = 'drop'
      )
    
    recurrence_status <- merged_sing_df %>%
      group_by(Pathway, Timepoint, study) %>%
      summarise(
        `Comparison Type` = "Recurrence Status",
        p_value = tryCatch(
          wilcox.test(Singscore ~ recurrence_status)$p.value,
          error = function(e) NA
        ),
        .groups = 'drop'
      )
    
    # Combine and sort
    bind_rows(response_status, recurrence_status) %>%
      arrange(p_value)
  })
  
  # Store the current pagination state
  current_display_count <- reactiveVal(10)
  
  # Update the number of rows to display when DataTable state changes
  observeEvent(input$globalTopSignificantTable_state, {
    state <- input$globalTopSignificantTable_state
    if (!is.null(state$length) && state$length > 0) {
      current_display_count(state$length)
    }
  })
  
  # Render Global Top Significant Comparisons
  output$globalTopSignificantTable <- DT::renderDataTable({
    req(global_comparisons())
    
    # Display the appropriate number of rows
    DT::datatable(global_comparisons()[1:current_display_count(), ],
                  options = list(pageLength = current_display_count()))
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
    
    if (input$timepoint != "Both") {
      data <- data %>% filter(Timepoint == input$timepoint)
    }
    
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
  
  output$sampleDistributionPlot <- renderPlot({
    req(selected_data())
    
    # Filter and prepare the data
    data <- selected_data()
    
    # Dynamically adjust the grouping based on comparison type
    if (input$comparison == "Response") {
      data <- data %>%
        group_by(study, Timepoint, Response) %>%
        summarise(Count = n(), .groups = 'drop')
      
      # Plot
      ggplot(data, aes(x = study, y = Count, fill = Response)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
        facet_wrap(~ Timepoint, scales = "free_y") +
        geom_text(aes(label = Count), 
                  position = position_dodge(width = 0.7), 
                  vjust = -0.5, 
                  size = 4) +
        theme_minimal() +
        labs(
          title = "Sample Distribution by Study, Timepoint, and Response",
          x = "Study",
          y = "Sample Count",
          fill = "Response"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.spacing = unit(0.5, "lines"),
          strip.background = element_rect(fill = "grey80", color = "black"),
          strip.text = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)
        )
    } else {
      data <- data %>%
        group_by(study, Timepoint, Comparison) %>%
        summarise(Count = n(), .groups = 'drop')
      
      # Plot
      ggplot(data, aes(x = study, y = Count, fill = Comparison)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
        facet_wrap(~ Timepoint, scales = "free_y") +
        geom_text(aes(label = Count), 
                  position = position_dodge(width = 0.7), 
                  vjust = -0.5, 
                  size = 4) +
        theme_minimal() +
        labs(
          title = "Sample Distribution by Study, Timepoint, and Recurrence Status",
          x = "Study",
          y = "Sample Count",
          fill = "Recurrence Status"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.spacing = unit(0.5, "lines"),
          strip.background = element_rect(fill = "grey80", color = "black"),
          strip.text = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)
        )
    }
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
    DT::datatable(selected_data()[, !(names(selected_data()) %in% c("MPRvNMPR", "X"))])
  })
  
  # Download all plots as a ZIP
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("Pathway_Enrichment_Plots.zip")
    },
    content = function(file) {
      # Create a temporary directory
      temp_dir <- tempdir()
      file_paths <- c()
      
      # Loop through each selected pathway and generate the plot
      lapply(input$pathway, function(path) {
        # Sanitize the filename: remove or replace problematic characters
        safe_path <- gsub("[:/\\*?\"<>|]", "-", path)  # replaces forbidden characters with "-"
        safe_path <- gsub("\\s+", "_", safe_path)      # replace spaces with underscores
        
        # Ensure it's unique
        plot_path <- file.path(temp_dir, paste0(safe_path, ".png"))
        
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
        
        # Generate the plot
        plot <- ggplot(data, aes(x = Timepoint, y = Singscore, fill = Comparison)) +
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
        
        # Save the plot as PNG
        ggsave(filename = plot_path, plot = plot, width = 10, height = 5, bg = "white", dpi = 300)
        
        # Add to the list of files
        file_paths <<- c(file_paths, plot_path)
      })
      
      # Zip all the files
      zip::zipr(zipfile = file, files = file_paths, recurse = FALSE)
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