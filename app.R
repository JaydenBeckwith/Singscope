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
library(shinythemes)
library(shinyjs)
library(plotly)
library(heatmaply)

merged_sing_df <- readRDS("merged_sing_df.rds")
gmt_data <- clusterProfiler::read.gmt("data/20241021_188genelist_withphenotypes.gmt")
gmt_data
# === Define UI ===
ui <- fluidPage(
  theme = shinytheme("sandstone"),
  useShinyjs(),
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
        )
      ),
      
      fluidRow(
        column(
          12,
          actionButton("showExample", "Show Example Data Format", icon = icon("eye")),
          div(
            id = "exampleData",
            style = "display: none;", # Hide initially
            br(),
            h4("Example Gene Expression Matrix"),
            DT::dataTableOutput("exampleExprMatrix"),
            br(),
            h4("Example Metadata"),
            DT::dataTableOutput("exampleMetadata")
          )
        )
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
                      choices = c("Response Status" = "Response", 
                                  "Recurrence Status" = "Recurrence Status",
                                  "Dynamics by Response" = "Dynamics_Response",
                                  "Dynamics by Recurrence" = "Dynamics_Recurrence"
                      )),
          selectInput("timepoint", "Select Timepoint", choices = c("Both", "Baseline", "Week 6")),
          fluidRow(
            column(12, 
                   h4("Create Custom Cohort(s)"),
                   textInput("cohortNameInput", "Cohort Name", value = "MyCohort")
            )
          ),
          fluidRow(
            column(12, 
                   h4("Select Existing Cohort"),
                   selectInput("selectedCohort", "Select Cohort", choices = NULL, multiple = FALSE)
            )
          ),
          fluidRow(
            column(12,
                   h5(textOutput("cohortCount"))
            )
          ),
          fluidRow(
            column(12,
                   div(
                     style = "display: flex;  gap: 10px; margin-bottom: 10px;",
                     actionButton("createCohort", "Create Cohort"),
                     actionButton("deleteCohort", "Delete Cohort")
                   )
            )
          ),
          fluidRow(
            column(12,
                   div(
                     style = "display: flex;  gap: 10px;  margin-bottom: 10px;",
                     downloadButton("downloadPlot", "Download Plots"),
                     downloadButton("downloadTable", "Download Table")
                   )
            )
          ),
          fluidRow(
            column(12,
                   div(
                     style = "display: flex;  gap: 10px;  margin-bottom: 10px;",
                     downloadButton("downloadSingscoreMatrix", "Download Singscore Matrix")
                   )
            )
          ),
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
                   h4("Select on rows to create a custom cohort."),
                   DT::dataTableOutput("dataTable")
            )
          ),
          fluidRow(
            h3("Selected Signature Plots"),
            column(12, uiOutput("dynamicPlots"))
          )
        )
      )
    ),
    
    # New Tab for Pathway Correlation Analysis
    tabPanel(
      "Signature Correlation Analysis",
      sidebarLayout(
        sidebarPanel(
          selectInput("correlationMethod", "Select Correlation Method",
                      choices = c("pearson", "spearman", "kendall")),
          selectInput("timepointFilter", "Select Timepoint", 
                      choices = c("Both", "Baseline", "Week 6")),
          selectInput("cohortFilter", "Select Cohort", 
                      choices = c("All", unique(merged_sing_df$study))),
          actionButton("computeCorrelation", "Compute Correlation"),
          br(),
          downloadButton("downloadCorrelation", "Download Correlation Matrix")
        ),
        mainPanel(
          fluidRow(
            column(12,
                   h3("Pathway Correlation Scatter Plot"),
                   plotOutput("pathwayCorrelationScatter"),
                   div(style = "margin-bottom: 20px;")
            )
          ),
          fluidRow(
            column(12,
                   h3("Pathway Correlation Heatmap"),
                   plotlyOutput("correlationHeatmap", height = "600px"),
                   div(style = "margin-bottom: 20px;")
            )
          )
        )
      )
    )))




# === Define Server Logic ===
server <- function(input, output, session) {
  
  
  # Example DataFrames
  example_expr <- data.frame(
    Gene = c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419"),
    CF27130 = c(991, 0, 524),
    CF27131 = c(205, 0, 582),
    CF27132 = c(266, 17, 383),
    CF27133 = c(446, 2, 570),
    CF27134 = c(434, 0, 482)
  )
  
  example_meta <- data.frame(
    sample_id = c("CF27130", "CF27131", "CF27132", "CF27133", "CF27134"),
    study = c("OpacinNeo", "OpacinNeo", "OpacinNeo", "OpacinNeo", "OpacinNeo"),
    recurrence_status = c(0, 0, 0, 0, 0),
    Gender = c("M", "F", "M", "F", "M"),
    MPRvNMPR = c(1, 0, 1, 1, 1),
    Timepoint = c("Baseline", "Baseline", "Baseline", "Baseline", "Baseline")
  )
  
  # Render example data
  output$exampleExprMatrix <- DT::renderDataTable({
    DT::datatable(example_expr, options = list(pageLength = 5))
  })
  
  output$exampleMetadata <- DT::renderDataTable({
    DT::datatable(example_meta, options = list(pageLength = 5))
  })
  
  # Toggle visibility on button click
  observeEvent(input$showExample, {
    toggle("exampleData")
  })
  

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
      
      # Update the global reference in memory
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
  
  
  # Reactive expression to get associated genes from GMTy
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
          input$comparison == "Recurrence Status" ~ ifelse(recurrence_status == 1, "Recurrence", "No Recurrence"),
          input$comparison == "Dynamics_Response" ~ ifelse(MPRvNMPR == 1, "MPRs", "NMPRs"),
          input$comparison == "Dynamics_Recurrence" ~ ifelse(recurrence_status == 1, "Recurrence", "No Recurrence")
        )
      )
    
    # If a cohort is selected, filter by it
    cohort_name <- input$selectedCohort
    if (!is.null(cohort_name) && cohort_name != "") {
      cohort_data <- cohorts$groups[[cohort_name]]
      if (!is.null(cohort_data)) {
        data <- data %>% filter(sample_id %in% cohort_data$sample_id)
      }
    }
    
    # For Dynamics options, we filter only Baseline and Week 6
    if (input$comparison %in% c("Dynamics_Response", "Dynamics_Recurrence")) {
      data <- data %>%
        filter(Timepoint %in% c("Baseline", "Week 6")) %>%
        group_by(Pathway, study, Comparison, Timepoint, sample_id) %>%
        distinct() %>%
        ungroup()
    }
    
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
    
    # === Deduplicate samples to avoid double counting ===
    data <- data %>%
      distinct(study, Timepoint, Comparison, sample_id, .keep_all = TRUE)
    
    # === SWITCH LOGIC FOR DYNAMICS OR STANDARD COMPARISON ===
    if (input$comparison %in% c("Response", "Dynamics_Response")) {
      data <- data %>%
        group_by(study, Timepoint, Response) %>%
        summarise(Count = n(), .groups = 'drop')
      
      # Plot for Response
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
      
    } else if (input$comparison %in% c("Recurrence Status", "Dynamics_Recurrence")) {
      data <- data %>%
        group_by(study, Timepoint, Comparison) %>%
        summarise(Count = n(), .groups = 'drop')
      
      # Plot for Recurrence
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
  
  
  observe({
    lapply(input$pathway, function(path) {
      output[[paste0("plot_", path)]] <- renderPlot({
        
        # Filter data for the selected pathway
        data <- selected_data() %>% filter(Pathway == path)
        
        # --- DEBUGGING POINT ---
        print(paste0("Pathway: ", path))
        print(head(data))
        print(input$comparison)
        
        # If there is no data, let's break early
        if (nrow(data) == 0) {
          showNotification("No data available for this pathway.", type = "error")
          return(NULL)
        }
        
        # --- SWITCH LOGIC FOR DYNAMICS VS STANDARD COMPARISON ---
        if (input$comparison %in% c("Dynamics_Response", "Dynamics_Recurrence")) {
          
          # Compute Wilcoxon p-values based on Timepoint for Dynamics
          p_values_df <- tryCatch({
            data %>%
              group_by(study, Comparison) %>%
              summarise(
                p_value = wilcox.test(Singscore ~ Timepoint)$p.value,
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
          }, error = function(e) {
            print(paste0("Wilcoxon Error: ", e$message))
            return(NULL)
          })
          
          if (is.null(p_values_df)) {
            showNotification("Wilcoxon test failed.", type = "error")
            return(NULL)
          }
          
          # --- PLOT FOR DYNAMICS ---
          ggplot(data, aes(x = Timepoint, y = Singscore, fill = Timepoint)) +
            geom_boxplot(width = 0.5, alpha = 0.8, position = position_dodge(0.7), outlier.shape = NA) +
            geom_jitter(
              shape = 21, size = 3,
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
              color = "black"
            ) +
            scale_fill_manual(values = c("Baseline" = "darkblue", "Week 6" = "orange")) +
            facet_grid(Comparison ~ study) +
            geom_text(
              data = p_values_df,
              aes(x = 1.5, y = max(data$Singscore, na.rm = TRUE) * 0.95, label = label),
              inherit.aes = FALSE,
              size = 5
            ) +
            theme_minimal() +
            labs(
              title = paste(gsub("_", " ", path), "Enrichment - Dynamics"),
              x = "Time",
              y = paste(path, "Singscore"),
              fill = "Time"
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
          
        } else {
          # --- STANDARD RESPONSE OR RECURRENCE PLOTTING ---
          
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
        }
      })
    })
  })
  
  # Render Data Table
  output$dataTable <- DT::renderDataTable({
  DT::datatable(
    selected_data()[, !(names(selected_data()) %in% c("MPRvNMPR", "X"))], 
    options = list(
      columnDefs = list(
        list(targets = 0, checkboxes = TRUE)  # Adds checkboxes to the first column
      )
    ),
    selection = list(
      mode = 'multiple',
      target = 'row'
    ),
    rownames = FALSE,
    class = 'display'
  )
})
  
  cohorts <- reactiveValues(groups = list())
  
  # === Create Cohort ===
  observeEvent(input$createCohort, {
    selected_rows <- input$dataTable_rows_selected
    
    if (length(selected_rows) == 0) {
      showNotification("No samples selected to create a cohort.", type = "error")
      return(NULL)
    }
    
    cohort_name <- input$cohortNameInput
    if (cohort_name %in% names(cohorts$groups)) {
      showNotification("Cohort name already exists. Choose a different name.", type = "error")
      return(NULL)
    }
    
    # Subset the selected samples from the DataTable
    cohort_data <- selected_data()[selected_rows, ]
    cohorts$groups[[cohort_name]] <- cohort_data
    
    # Update dropdown
    updateSelectInput(session, "selectedCohort", choices = names(cohorts$groups))
    showNotification(paste("Cohort", cohort_name, "created successfully!"), type = "message")
  })
  
  # === Delete Cohort ===
  observeEvent(input$deleteCohort, {
    cohort_name <- input$selectedCohort
    if (cohort_name %in% names(cohorts$groups)) {
      cohorts$groups[[cohort_name]] <- NULL
      updateSelectInput(session, "selectedCohort", choices = names(cohorts$groups))
      showNotification(paste("Cohort", cohort_name, "deleted successfully!"), type = "message")
    } else {
      showNotification("Cohort not found.", type = "error")
    }
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
        
        # Skip if data is empty
        if (nrow(data) == 0) {
          showNotification(paste("No data for pathway:", path), type = "error")
          return(NULL)
        }
        
        # === SWITCH LOGIC FOR DYNAMICS OR STANDARD COMPARISON ===
        if (input$comparison %in% c("Dynamics_Response", "Dynamics_Recurrence")) {
          
          # Compute Wilcoxon p-values based on Timepoint for Dynamics
          p_values_df <- tryCatch({
            data %>%
              group_by(study, Comparison) %>%
              summarise(
                p_value = wilcox.test(Singscore ~ Timepoint)$p.value,
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
          }, error = function(e) {
            print(paste("Wilcoxon Error:", e$message))
            return(NULL)
          })
          
          # Generate the plot for Dynamics
          plot <- ggplot(data, aes(x = Timepoint, y = Singscore, fill = Timepoint)) +
            geom_boxplot(width = 0.5, alpha = 0.8, position = position_dodge(0.7), outlier.shape = NA) +
            geom_jitter(
              shape = 21, size = 3,
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
              color = "black"
            ) +
            scale_fill_manual(values = c("Baseline" = "darkblue", "Week 6" = "orange")) +
            facet_grid(Comparison ~ study) +
            geom_text(
              data = p_values_df,
              aes(x = 1.5, y = max(data$Singscore, na.rm = TRUE) * 0.95, label = label),
              inherit.aes = FALSE,
              size = 5
            ) +
            theme_minimal() +
            labs(
              title = paste(gsub("_", " ", path), "Enrichment - Dynamics"),
              x = "Time",
              y = paste(path, "Singscore"),
              fill = "Time"
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
          
        } else {
          # === STANDARD RESPONSE OR RECURRENCE PLOTTING ===
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
        }
        
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
  
  output$downloadSingscoreMatrix <- downloadHandler(
    filename = function() {
      cohort_name <- ifelse(is.null(input$selectedCohort) || input$selectedCohort == "", "All_Cohorts", input$selectedCohort)
      timepoint <- input$timepoint
      paste0(cohort_name, "_", timepoint, "_Singscore_Matrix.csv")
    },
    content = function(file) {
      
      # === 1 Get the Data
      cohort_name <- gsub(" \\(.*\\)$", "", input$selectedCohort)
      
      if (!is.null(cohort_name) && cohort_name %in% names(cohorts$groups)) {
        cohort_data <- cohorts$groups[[cohort_name]]
      } else {
        cohort_data <- selected_data()
      }
      
      # ===  Filter by Timepoint
      if (input$timepoint != "Both") {
        cohort_data <- cohort_data %>% filter(Timepoint == input$timepoint)
      }
      
      # ===  Create the Singscore Matrix
      singscore_matrix <- cohort_data %>%
        select(sample_id, Pathway, Singscore) %>%
        pivot_wider(names_from = Pathway, values_from = Singscore) %>%
        arrange(sample_id) %>%
        column_to_rownames(var = "sample_id")
      
      # ===  Save as CSV
      write.csv(singscore_matrix, file, row.names = TRUE)
    }
  )
  
  observeEvent(input$computeCorrelation, {
    req(input$cohortFilter, input$timepointFilter)
    
    showNotification("Computing pathway correlations...", type = "message")
    
    # === Get the appropriate cohort
    cohort_name <- input$cohortFilter
    if (cohort_name == "All") {
      cohort_data <- merged_sing_df
    } else {
      cohort_data <- cohorts$groups[[cohort_name]]
    }
    
    # === Filter by Timepoint
    if (input$timepointFilter != "Both") {
      cohort_data <- cohort_data %>% filter(Timepoint == input$timepointFilter)
    }
    
    # === Create a Singscore Matrix
    singscore_matrix <- cohort_data %>%
      select(sample_id, Pathway, Singscore) %>%
      pivot_wider(names_from = Pathway, values_from = Singscore) %>%
      column_to_rownames("sample_id")
    
    # === Enforce Valid Method and Clean Whitespace ===
    cor_method <- trimws(input$correlationMethod)  # Remove any extra spaces
    cor_method <- match.arg(cor_method, c("pearson", "kendall", "spearman"))
    
    # === Compute the Correlation Matrix
    correlation_matrix <- cor(singscore_matrix, method = cor_method, use = "pairwise.complete.obs")
    
    showNotification("Rendering interactive heatmap...", type = "message")
    
    # === Plot the Interactive Heatmap
    output$correlationHeatmap <- renderPlotly({
      heatmaply::heatmaply(
        correlation_matrix,
        main = "Interactive Pathway Correlation Heatmap",
        xlab = "Pathways",
        ylab = "Pathways",
        scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
          low = "blue", mid = "white", high = "red", midpoint = 0
        ),
        dendrogram = "both",
        k_row = 3, k_col = 3,
        grid_color = "grey50",
        grid_width = 0.2,
        label_names = c("Pathway 1", "Pathway 2", "Correlation"),
        fontsize_row = 8,
        fontsize_col = 8,
        showticklabels = c(TRUE, TRUE)
      )
    })
    
    showNotification("Interactive heatmap rendering complete!", type = "message")
    
    # === Allow Download of the Correlation Matrix
    output$downloadCorrelation <- downloadHandler(
      filename = function() {
        paste0("Pathway_Correlation_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(correlation_matrix, file, row.names = TRUE)
      }
    )
  })
  
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