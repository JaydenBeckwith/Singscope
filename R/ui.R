library(shiny)
library(shinythemes)
library(shinyjs)
library(shinyWidgets)
library(plotly)
library(DT)
library(shinycssloaders)


# === Load Data ===
rds_path <- if (file.exists("merged_sing_df.rds")) {
  "merged_sing_df.rds"
} else {
  "/srv/shiny-server/data/merged_sing_df.rds"
}
merged_sing_df <- readRDS(rds_path)

ui <- navbarPage(
   div("Singscope", style = "font-size: 24px; font-weight: bold;"),
  theme = shinytheme("sandstone"),
  id = "mainTabs",
  useShinyjs(),
  tags$head(
    tags$style(HTML(".hamburger-btn {
      background-color: #F5F5F5;
      border: 2px solid #c2b49a;
      color: #c2b49a;
      font-size: 20px;
      padding: 6px 12px;
      margin: 15px 0 10px 10px;
      border-radius: 4px;
    }
    .hamburger-btn:hover {
      background-color: #fff;
      color: #a08760;
      border-color: #a08760;
    }")),
    ### starts the zoom at 80% upon boot up 
    tags$script(HTML("
      document.addEventListener('DOMContentLoaded', function() {
        document.body.style.zoom = '80%';
      });
    "))
  ),
  
  tabPanel(
    "Data Import",
    fluidRow(
      column(
        6, offset = 3,
        wellPanel(
          h3("Import Data"),
          fileInput("exprMatrix", "Upload Gene Expression Matrix (.csv, .tsv)", accept = c(".csv", ".tsv")),
          fileInput("singMatrix", "Upload Singscore Matrix (.csv, .tsv)", accept = c(".csv", ".tsv")),
          fileInput("metadata", "Upload Clinical Data (.csv, .tsv)", accept = c(".csv", ".tsv")),
          checkboxInput("mergeToExample", "Merge with example dataset", value = TRUE),
          actionButton("submitData", "Submit Data"),
          br(),
          verbatimTextOutput("preprocessLog"),
          div(style = "margin-bottom: 20px;")
        )
      )
    ),
    fluidRow(
    column(
      12,
      div(
        style = "display: flex; gap: 10px;",
        actionButton("showExample", "Show Example Data Format", icon = icon("eye")),
        actionButton("useExampleData", "Use Example Data", icon = icon("upload"))
      ),
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
      column(12, actionButton("toggleSignatureSidebar", HTML("&#9776;"), class = "hamburger-btn"))
    ),
    sidebarLayout(
      div(id = "signatureSidebar",
          sidebarPanel(
            textOutput("signatureCountText"), br(),
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
        h3("Top Statistically Significant Comparisons"),
        DT::dataTableOutput("globalTopSignificantTable"),
        br(),
        fluidRow(
          column(12,
                h3("Sample Distribution"),
                withSpinner(plotOutput("sampleDistributionPlot", height = "500px"))
          )
        ),
        fluidRow(
          column(12,
                h3("Mutation Distribution"),
                withSpinner(plotlyOutput("mutationPieChart", height = "500px"))
          )
        ),
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
      column(12, actionButton("toggleCorrelationSidebar", HTML("&#9776;"), class = "hamburger-btn"))
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
            tags$hr(),
            downloadButton("downloadCorrelation", "Download Correlation Matrix")
          )
      ),
      mainPanel(
        fluidRow(
          column(12,
                 h3("Correlation Signatures"),
                 withSpinner(plotlyOutput("correlationHeatmap", height = "1000px", width = "100%"))
          )
        ),
        fluidRow(
          column(12,
                 h3("Change in Correlation Across Time"),
                 withSpinner(plotlyOutput("deltaCorrelationHeatmap", height = "500px", width = "100%"))
          )
        )
      )
    )
  ),
  
  tabPanel(
    "Survival Analysis",
    fluidRow(
      column(12, actionButton("toggleSurvivalSidebar", HTML("&#9776;"), class = "hamburger-btn"))
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