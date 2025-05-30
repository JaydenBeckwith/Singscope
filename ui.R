
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
ui <- ui <- fluidPage(
  theme = shinytheme("sandstone"),
  useShinyjs(),
  tags$head(
    tags$style(HTML(".hamburger-btn {
      background-color: #F5F5F5; /* sandstone-like light tone */
      border: 2px solid #c2b49a; /* sandstone accent */
      color: #c2b49a;
      font-size: 20px;
      padding: 6px 12px;
      margin: 15px 0 10px 10px;  /* space around button */
      border-radius: 4px;
    }
    .hamburger-btn:hover {
      background-color: #fff;
      color: #a08760;
      border-color: #a08760;
    }"))
  ),
  titlePanel("Signature Enrichment Comparison"),
  
  tabsetPanel(
    tabPanel(
      "Data Import",
      fluidRow(
        column(
          6, offset = 3,
          wellPanel(
            h3("Import Data"),
            fileInput("exprMatrix", "Upload Gene Expression Matrix (.csv, .tsv)", accept = c(".csv", ".tsv")),
            fileInput("metadata", "Upload Metadata (.csv, .tsv)", accept = c(".csv", ".tsv")),
            actionButton("submitData", "Submit Data"),
            div(style = "margin-bottom: 20px;")
          )
        )
      ),
      fluidRow(
        column(
          12,
          actionButton("showExample", "Show Example Data Format", icon = icon("eye")),
          div(id = "exampleData", style = "display: none;",
              br(), h4("Example Gene Expression Matrix"), DT::dataTableOutput("exampleExprMatrix"),
              br(), h4("Example Metadata"), DT::dataTableOutput("exampleMetadata"))
        )
      ),
      fluidRow(
        column(12,
               h3("Uploaded Data Preview"),
               DT::dataTableOutput("previewExprMatrix"),
               DT::dataTableOutput("previewMetadata")
        )
      )
    ),
    
    tabPanel(
      "Signature Analysis",
      fluidRow(
        column(12,
               actionButton("toggleSignatureSidebar", HTML("&#9776;"), class = "hamburger-btn")
        )
      ),
      sidebarLayout(
        div(id = "signatureSidebar",
            sidebarPanel(
              pickerInput("pathway", "Select Signatures", choices = unique(merged_sing_df$Pathway),
                          selected = unique(merged_sing_df$Pathway)[1], multiple = TRUE, options = list(`actions-box` = TRUE)),
              selectInput("study", "Select Study", choices = c("All", unique(merged_sing_df$study))),
              selectInput("comparison", "Select Comparison Type", choices = c("Response Status" = "Response", "Recurrence Status" = "Recurrence Status", "Dynamics by Response" = "Dynamics_Response", "Dynamics by Recurrence" = "Dynamics_Recurrence")),
              selectInput("timepoint", "Select Timepoint", choices = c("Both", "Baseline", "Week 6")),
              textInput("cohortNameInput", "Cohort Name", value = "MyCohort"),
              selectInput("selectedCohort", "Select Cohort", choices = NULL, multiple = FALSE),
              textOutput("cohortCount"),
              div(style = "display: flex; gap: 10px; margin-bottom: 10px;",
                  actionButton("createCohort", "Create Cohort"),
                  actionButton("deleteCohort", "Delete Cohort")),
              div(style = "display: flex; gap: 10px; margin-bottom: 10px;",
                  downloadButton("downloadPlot", "Download Plots"),
                  downloadButton("downloadTable", "Download Table")),
              div(style = "margin-bottom: 10px;",
                  downloadButton("downloadSingscoreMatrix", "Download Singscore Matrix")),
              h4("Genes in Selected Signature(s)"),
              uiOutput("geneList"),
              width = 3
            )
        ),
        mainPanel(
          h3("Global Top Statistically Significant Comparisons"),
          DT::dataTableOutput("globalTopSignificantTable"),
          br(),
          fluidRow(
            column(7, h3("Sample Distribution"), plotOutput("sampleDistributionPlot", height = "500px")),
            column(5, h3("Mutation Distribution"), plotlyOutput("mutationPieChart", height = "600px"))
          ),
          br(),
          h3("Selected Sample Data"),
          h4("Select on rows to create a custom cohort."),
          DT::dataTableOutput("dataTable"),
          h3("Selected Signature Plots"),
          uiOutput("dynamicPlots")
        )
      )
    ),
    
    tabPanel(
      "Signature Correlation Analysis",
      fluidRow(
        column(12,
               actionButton("toggleCorrelationSidebar", HTML("&#9776;"), class = "hamburger-btn")
        )
      ),
      sidebarLayout(
        div(id = "correlationSidebar",
            sidebarPanel(
              h4("Correlation Settings"),
              selectInput("correlationMethod", "Correlation Method", choices = c("pearson", "spearman", "kendall")),
              selectInput("timepointFilter", "Timepoint", choices = c("Both", "Baseline", "Week 6")),
              selectInput("cohortFilter", "Cohort", choices = c("All", unique(merged_sing_df$study))),
              actionButton("computeCorrelation", "Compute Correlation", class = "btn-primary"),
              tags$hr(), h4("Temporal Trajectory Analysis"),
              pickerInput("trajectoryPathways", "Select Pathways", choices = unique(merged_sing_df$Pathway), multiple = TRUE, options = list(`actions-box` = TRUE)),
              selectInput("studyFilter", "Study", choices = c("All", unique(merged_sing_df$study))),
              actionButton("computeTrajectory", "Compute Trajectory", class = "btn-info"),
              tags$hr(), downloadButton("downloadCorrelation", "Download Correlation Matrix")
            )
        ),
        mainPanel(
          fluidRow(
            column(12,
                   h3("Correlation Signatures"),
                   plotlyOutput("correlationHeatmap", height = "1000px", width = "100%")
            )
          ),
          fluidRow(
            column(12,
                   h3("Change in Correlation Across Time"),
                   plotlyOutput("deltaCorrelationHeatmap", height = "500px", width = "100%")
            )
          )
        )
      )
    ),
    
    tabPanel(
      "Survival Analysis",
      fluidRow(
        column(12,
               actionButton("toggleSurvivalSidebar", HTML("&#9776;"), class = "hamburger-btn")
        )
      ),
      sidebarLayout(
        div(id = "survivalSidebar",
            sidebarPanel(
              selectInput("survivalTime", "Select Time Variable", choices = c("time_to_event")),
              selectInput("survivalEvent", "Select Event Variable", choices = c("event")),
              selectInput("groupingVariable", "Group By", choices = c("Mutation", "Response", "Custom Group")),
              selectInput("timepointSurv", "Filter by Timepoint", choices = c("Both", "Baseline", "Week 6")),
              selectInput("studySurv", "Filter by Study", choices = c("All", unique(merged_sing_df$study))),
              actionButton("runSurvival", "Run Survival Analysis")
            )
        ),
        mainPanel(
          h3("Kaplan-Meier Survival Curve"),
          plotlyOutput("kmPlot", height = "500px")
        )
      )
    )
  )
)

