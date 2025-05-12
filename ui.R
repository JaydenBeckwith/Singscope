
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
          pickerInput("trajectoryPathways", "Select Pathways for Trajectory",
                      choices = unique(merged_sing_df$Pathway),
                      multiple = TRUE,
                      options = list(`actions-box` = TRUE)),
          selectInput("studyFilter", "Select Study", 
                      choices = c("All", unique(merged_sing_df$study))),
          actionButton("computeTrajectory", "Compute Trajectory"),
          br(),
          downloadButton("downloadCorrelation", "Download Correlation Matrix")
        ),
        mainPanel(
          fluidRow(
            column(12,
                   h3("Pathway Correlation Heatmap"),
                   plotlyOutput("correlationHeatmap", height = "1000px", width = "1000px"),
                   div(style = "margin-bottom: 20px;")
            )
          ),
          mainPanel(
            fluidRow(
              column(
                12,  # 50% width of the row
                h3("Pathway Correlation Trajectories"),
                plotlyOutput("trajectoryPlot", height = "500px")
              ),
              column(
                12,  # 50% width of the row
                h3("Delta Correlation Heatmap"),
                plotlyOutput("deltaCorrelationHeatmap", height = "500px")
              )
            ),
            div(style = "margin-bottom: 20px;")
          )
        )
      )
    )))
