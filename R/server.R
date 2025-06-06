library(GSEABase)
source("preprocess_data.R")

# === Define Server Logic ===
EXAMPLE_DATA <- if (file.exists("example_data_neotrio_neopele.rds")) {
  "merged_sing_df.rds"
} else {
  "/srv/shiny-server/data/example_data_neotrio_neopele.rds"
}
merged_sing_df <- readRDS(EXAMPLE_DATA)


gmt_path <- if (file.exists("/srv/shiny-server/data/20251505_240genelist_withphenotypes.gmt")) {
  "/srv/shiny-server/data/20251505_240genelist_withphenotypes.gmt"
} else {
  "/srv/shiny-server/data/20251505_240genelist_withphenotypes.gmt"
}

gmt_data_show <- read.csv("/srv/shiny-server/data/singscore_gene_enrichment_list.csv")

server <- function(input, output, session) {

  version <- reactive({
    readLines("www/version.txt", warn = FALSE)[1]
  })

  last_updated <- reactive({
    file.info("www/version.txt")$mtime |> format("%d %b %Y %H:%M")
  })

  output$app_version <- renderUI({
    HTML(glue::glue('
      <a href="www/version.txt"
         title="Last updated: {last_updated()}"
         target="_blank"
         style="
           text-decoration: none;
           background-color: #e6f3f7;
           color: #29A3C1;
           font-size: 12px;
           padding: 4px 10px;
           border-radius: 4px;
           font-weight: bold;
           margin-right: 10px;
           border: 1px solid #29A3C1;
           box-shadow: 1px 1px 3px rgba(0,0,0,0.1);
         ">
         {version()}
      </a>
    '))
  })

  tryCatch({
  # Example DataFrames
  example_expr <- data.frame(
    ensembl_gene_id = c("ENSG00000000003", "ENSG00000000005", "ENSG00000000419"),
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
    Mutation = c("BRAFV600E", "BRAFV600E", "BRAFV600E", "BRAFV600E", "BRAFV600K"),
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
  
  expr_data <- reactiveVal(NULL)
  metadata_data <- reactiveVal(NULL)
  
  observeEvent(input$submitData, {
    req(input$exprMatrix, input$metadata)

    shinyjs::disable("submitData")
    shinyjs::html("submitData", '<i class="fa fa-spinner fa-spin"></i> Loading...')

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

    log_text <- reactiveVal("")

    append_log <- function(msg) {
      isolate(log_text(paste(log_text(), msg, sep = "\n")))
    }

    output$preprocessLog <- renderText({
      log_text()
    })

    # === 4. Preprocessing Status Update ===
    output$preprocessStatus <- renderUI({
      div(
        h4("Preprocessing in progress..."),
        p("Please wait while your data is being processed.")
      )
    })

    tryCatch({
      log_text("Starting preprocessing...")  

      new_data <- preprocess_data(
        exprMatrixPath = input$exprMatrix$datapath,
        metadataPath = input$metadata$datapath,
        cohortName = input$cohortName,
        gmtPath = gmt_path,
        existing_data_merge = input$mergeToExample,
        example_data_path = EXAMPLE_DATA,
        logger = append_log
      )

      merged_sing_df <<- new_data

      updatePickerInput(session, "pathway", choices = unique(merged_sing_df$Pathway))
      updateSelectInput(session, "study", choices = c("All", unique(merged_sing_df$study)))

      showNotification("Data imported and merged successfully!", type = "message")
      append_log("Preprocessing complete.")

    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      append_log(paste("Error:", e$message))

    }, finally = {
      shinyjs::enable("submitData")
      shinyjs::html("submitData", "Submit Data")
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
      genes <- gmt_data_show$gene[gmt_data_show$term == pathway]
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
        Comparison = dplyr::case_when(
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
    dplyr::bind_rows(response_status, recurrence_status) %>%
      arrange(p_value)
  })
  
  # Store the current pagination state
  #current_display_count <- reactiveVal(10)
  
  output$globalTopSignificantTable <- DT::renderDataTable({
  shiny::req(global_comparisons())
  
  DT::datatable(
    global_comparisons(),
    options = list(
      pageLength = 10,  # Default display
      lengthMenu = list(c(10, 25, 50, 100, -1), c('10', '25', '50', '100', 'All')),
      pagingType = "full_numbers"  # Adds first/last page buttons
    )
  )
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
    shiny::req(selected_data(), input$comparison)

    data <- selected_data()

    data <- data %>%
      distinct(study, Timepoint, Comparison, sample_id, .keep_all = TRUE)

    print(data)

    if (input$comparison %in% c("Response", "Dynamics_Response")) {
      validate(need("Response" %in% colnames(data), "Missing 'Response' column"))

    data <- data %>%
        group_by(study, Timepoint, Response) %>%
        summarise(Count = dplyr::n(), .groups = 'drop')
    message("error line?")
    #print(data)
      

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
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.spacing = unit(0.5, "lines"),
          strip.background = element_rect(fill = "grey80", color = "black"),
          strip.text = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
        )

    } else if (input$comparison %in% c("Recurrence Status", "Dynamics_Recurrence")) {
      data <- data %>%
        group_by(study, Timepoint, Comparison) %>%
        summarise(Count = dplyr::n(), .groups = 'drop')

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
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.spacing = unit(0.5, "lines"),
          strip.background = element_rect(fill = "grey80", color = "black"),
          strip.text = element_text(face = "bold", size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
        )
    }
  })

  #### mutation pie chart 
  output$mutationPieChart <- renderPlotly({
    shiny::req(selected_data())

    data <- selected_data()

    validate(
      need(all(c("sample_id", "Mutation") %in% colnames(data)), "Missing required columns for mutation pie chart")
    )

    mutation_summary <- data %>%
      mutate(PatientID = sub("_S[0-9]+$", "", sample_id)) %>%
      distinct(PatientID, Mutation) %>%
      group_by(Mutation) %>%
      summarise(Count = dplyr::n(), .groups = "drop")

    n_colors <- nrow(mutation_summary)
    palette <- RColorBrewer::brewer.pal(min(n_colors, 8), "Set2")
    if (n_colors > 8) {
      palette <- grDevices::colorRampPalette(palette)(n_colors)
    }

    plotly::plot_ly(
      mutation_summary,
      labels = ~Mutation,
      values = ~Count,
      type = 'pie',
      textinfo = 'label+percent',
      insidetextorientation = 'radial',
      marker = list(colors = palette)
    ) %>%
      layout(
        title = "",
        showlegend = TRUE,
        margin = list(l = 0, r = 0, b = 0, t = 30)
      )
  })

  #### nodal site pie chart 
  output$nodalSitePieChart <- renderPlotly({
    shiny::req(selected_data())
    
    data <- selected_data()
    
    validate(
      need(all(c("sample_id", "nodal_site") %in% colnames(data)), "Missing required columns for nodal site pie chart")
    )
    
    # 🔁 Map nodal site numbers to names
    nodal_map <- c("1" = "neck", "2" = "axilla", "3" = "groin")
    
    nodal_summary <- data %>%
      mutate(
        PatientID = sub("_S[0-9]+$", "", sample_id),
        nodal_site = as.character(nodal_site),
        nodal_site = nodal_map[nodal_site]
      ) %>%
      distinct(PatientID, nodal_site) %>%
      group_by(nodal_site) %>%
      summarise(Count = dplyr::n(), .groups = "drop")
    
    n_colors <- nrow(nodal_summary)
    palette <- RColorBrewer::brewer.pal(min(n_colors, 8), "Set2")
    if (n_colors > 8) {
      palette <- grDevices::colorRampPalette(palette)(n_colors)
    }
    
    plotly::plot_ly(
      nodal_summary,
      labels = ~nodal_site,
      values = ~Count,
      type = 'pie',
      textinfo = 'label+percent',
      insidetextorientation = 'radial',
      marker = list(colors = palette)
    ) %>%
      layout(
        title = "",
        showlegend = TRUE,
        margin = list(l = 0, r = 0, b = 0, t = 30)
      )
  })
  
  
  observe({
    lapply(input$pathway, function(path) {
      output[[paste0("plot_", path)]] <- renderPlot({
        
        # Filter data for the selected pathway
        data <- selected_data() %>% filter(Pathway == path)
        
        # --- DEBUGGING POINT ---
        print(paste0("Pathway: ", path))

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
                label = dplyr::case_when(
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
              axis.title.y = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.text.x = element_text(size = 14),  # adjust the number as needed
              axis.text.y = element_text(size = 14),  # optional: increase y-axis too
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              panel.spacing = unit(0.5, "lines"),
              strip.background = element_rect(fill = "grey80", color = "black"),
              strip.text = element_text(face = "bold", size = 14),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 14)
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
              label = dplyr::case_when(
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
              axis.title.y = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.text.x = element_text(size = 14),  # adjust the number as needed
              axis.text.y = element_text(size = 14),  # optional: increase y-axis too
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              panel.spacing = unit(0.5, "lines"),
              strip.background = element_rect(fill = "grey80", color = "black"),
              strip.text = element_text(face = "bold", size = 14),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 14)
            )
        }
      })
    })
  })
  
  # Render Selection Singscore Data Table
  output$dataTable <- DT::renderDataTable({
  dat <- selected_data()[, !(names(selected_data()) %in% c("MPRvNMPR", "X"))]
  dat <- cbind(" " = "", dat)  # Add blank first column for checkboxes

  DT::datatable(
    dat,
    extensions = "Select",
    selection = "none",
    options = list(
      pageLength = 10,  # Default display
      lengthMenu = list(c(10, 25, 50, 100, -1), c('10', '25', '50', '100', 'All')),
      pagingType = "full_numbers", 
      scrollX = TRUE, # Adds first/last page buttons
      columnDefs = list(
        list(
          targets = 0,
          orderable = FALSE,
          className = "select-checkbox"
        )
      ),
      select = list(
        style = "multi",
        selector = "td:first-child"
      ),
      dom = 'Blfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    ),
    rownames = FALSE,
    class = 'display',
    callback = JS(
      "
      // Add select all checkbox to header
      table.on('draw', function() {
        var checkbox = '<input type=\"checkbox\" id=\"select-all\">';
        if (!$('#select-all').length) {
          $(table.column(0).header()).html(checkbox);
        }
      });

      // Handle select all click
      $(document).on('click', '#select-all', function() {
        if (this.checked) {
          table.rows().select();
        } else {
          table.rows().deselect();
        }
      });

      // Keep Shiny in sync
      table.on('select deselect', function() {
        var indexes = table.rows({ selected: true }).indexes().toArray();
        Shiny.setInputValue('selected_rows', indexes);
      });
      "
    )
  )
})

cohorts <- reactiveValues(groups = list())

# === Create Cohort from Selected Rows ===
observeEvent(input$createCohort, {
  selected <- input$selected_rows
  
  if (is.null(selected) || length(selected) == 0) {
    showNotification("No samples selected to create a cohort.", type = "error")
    return(NULL)
  }

  cohort_name <- input$cohortNameInput
  if (cohort_name %in% names(cohorts$groups)) {
    showNotification("Cohort name already exists. Choose a different name.", type = "error")
    return(NULL)
  }

  # Subset data based on selected row indices
  cohort_data <- selected_data()[selected + 1, ]  # Add 1 since JS is 0-indexed
  cohorts$groups[[cohort_name]] <- cohort_data

  # Update dropdown menu
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
                label = dplyr::case_when(
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
              axis.title.y = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              panel.spacing = unit(0.5, "lines"),
              strip.background = element_rect(fill = "grey80", color = "black"),
              strip.text = element_text(face = "bold", size = 14),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 14)
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
              label = dplyr::case_when(
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
              axis.title.y = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              panel.spacing = unit(0.5, "lines"),
              strip.background = element_rect(fill = "grey80", color = "black"),
              strip.text = element_text(face = "bold", size = 14),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 14)
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
        tidyr::pivot_wider(names_from = Pathway, values_from = Singscore) %>%
        arrange(sample_id) %>%
        tibble::column_to_rownames(var = "sample_id")
      
      # ===  Save as CSV
      write.csv(singscore_matrix, file, row.names = TRUE)
    }
  )
  
  
  observeEvent(input$computeCorrelation, {
    req(input$cohortFilter, input$timepointFilter)
    
    showNotification("Computing pathway correlations...", type = "message")
    
    # === Get the appropriate cohort
    cohort_name <- input$cohortFilter
    timepoint_name <- input$timepointFilter
    
    if (cohort_name == "All") {
      cohort_data <- merged_sing_df
    } else {
      cohort_data <- merged_sing_df %>%
        filter(study == cohort_name)
    }
    
    # === Filter by Timepoint
    if (timepoint_name != "Both") {
      cohort_data <- cohort_data %>% filter(Timepoint == timepoint_name)
    }
    
    # === Create a Singscore Matrix
    singscore_matrix <- cohort_data %>%
      dplyr::select(sample_id, Pathway, Singscore) %>%
      tidyr::pivot_wider(names_from = Pathway, values_from = Singscore) %>%
      tibble::column_to_rownames("sample_id")
    
    # === Enforce Valid Method and Clean Whitespace ===
    cor_method <- trimws(input$correlationMethod)  
    cor_method <- match.arg(cor_method, c("pearson", "kendall", "spearman"))
    
    # === Compute the Correlation Matrix
    correlation_matrix <- cor(singscore_matrix, method = cor_method, use = "pairwise.complete.obs")
    
    showNotification("Rendering interactive heatmap...", type = "message")
    
    # === Plot the Interactive Heatmap
    output$correlationHeatmap <- renderPlotly({
      heatmaply::heatmaply(
        correlation_matrix,
        main = paste("Signature Correlation Heatmap -", cohort_name, "-", timepoint_name),
        xlab = NULL,
        ylab = NULL,
        scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
          low = "blue", mid = "white", high = "red", midpoint = 0
        ),
        dendrogram = "both",
        k_row = 3, k_col = 3,
        grid_color = NA,
        grid_width = 0.2,
        label_names = c("Pathway 1", "Pathway 2", "Correlation"),
        fontsize_row = 8,
        fontsize_col = 8,
        showticklabels = c(FALSE, FALSE)
      ) %>%
        plotly::layout(
          xaxis = list(
            showticklabels = FALSE,
            title = NULL,
            showline = FALSE,
            tickvals = list(),
            ticktext = list()
          ),
          yaxis = list(
            showticklabels = FALSE,
            title = NULL,
            showline = FALSE,
            tickvals = list(),
            ticktext = list()
          )
        )
    })
    
    showNotification("Interactive heatmap rendering complete!", type = "message")
    
    # === Allow Download of the Correlation Matrix
    output$downloadCorrelation <- downloadHandler(
      filename = function() {
        paste0("Pathway_Correlation_", cohort_name, "_", timepoint_name, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(correlation_matrix, file, row.names = TRUE)
      }
    )
  })
  
  observeEvent(input$cohortFilter, {
    if (input$cohortFilter != "") {
      shinyjs::click("computeCorrelation")
    }
  })
  
  observeEvent(input$computeTrajectory, {
    req(input$trajectoryPathways, input$studyFilter)
    
    showNotification("Computing trajectory of pathway correlations...", type = "message")
    
    # === Filter by study if not 'All'
    data <- if (input$studyFilter != "All") {
      merged_sing_df %>%
        filter(study == input$studyFilter)
    } else {
      merged_sing_df
    }
    
    # === Subset only the selected pathways
    data <- data %>%
      filter(Pathway %in% input$trajectoryPathways)
    
    # === Pivot to wide format to prepare for correlation calculation
    wide_data <- data %>%
      dplyr::select(sample_id, Pathway, Singscore, Timepoint) %>%
      tidyr::pivot_wider(names_from = Pathway, values_from = Singscore)
    
    # === Ensure numeric columns
    singscore_matrix <- wide_data %>%
      dplyr::select(-sample_id, -Timepoint) %>%
      mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
      as.data.frame()
    
    rownames(singscore_matrix) <- wide_data$sample_id
    timepoints_vector <- wide_data$Timepoint
    
    # === Compute correlation matrix per timepoint
    correlation_list <- list()
    for (tp in unique(timepoints_vector)) {
      idx <- which(timepoints_vector == tp)
      sub_matrix <- singscore_matrix[idx, , drop = FALSE]
      
      if (nrow(sub_matrix) > 1) {
        correlation_matrix <- cor(sub_matrix, use = "pairwise.complete.obs")
        correlation_list[[tp]] <- correlation_matrix
      }
    }
    
    # === Convert to long format
    correlation_df <- do.call(rbind, lapply(names(correlation_list), function(tp) {
      corr_mat <- correlation_list[[tp]]
      if (!is.null(corr_mat)) {
        df <- as.data.frame(as.table(corr_mat))
        df$Timepoint <- tp
        return(df)
      }
    }))
    colnames(correlation_df) <- c("Pathway1", "Pathway2", "Correlation", "Timepoint")
    
    selected_pairs <- combn(unique(input$trajectoryPathways), 2, simplify = FALSE)
    selected_pair_labels <- sapply(selected_pairs, function(x) paste(sort(x), collapse = " vs "))
    
    # Filter only selected pairs (avoid self-correlations)
    correlation_df <- correlation_df %>%
      filter(Pathway1 != Pathway2) %>%
      mutate(Pair = paste(pmin(Pathway1, Pathway2), pmax(Pathway1, Pathway2), sep = " vs ")) %>%
      filter(Pair %in% selected_pair_labels) %>%
      distinct(Pair, Timepoint, .keep_all = TRUE)
    
    # === Compute Delta Correlation Matrix (Week 6 - Baseline)
    if ("Baseline" %in% names(correlation_list) && "Week 6" %in% names(correlation_list)) {
      baseline_corr <- correlation_list[["Baseline"]]
      week6_corr <- correlation_list[["Week 6"]]
      delta_matrix <- week6_corr - baseline_corr
      
      # === Render Delta Correlation Heatmap
      output$deltaCorrelationHeatmap <- renderPlotly({
        heatmaply::heatmaply(
          delta_matrix,
          main = "Delta Correlation (Week 6 - Baseline)",
          xlab = "Pathways",
          ylab = "Pathways",
          colors = colorRampPalette(c("blue", "white", "red"))(100),
          dendrogram = "both"
        )
      })
    } else {
      showNotification("One of the timepoints (Baseline or Week 6) is missing.", type = "error")
    }
    
    showNotification("Trajectory analysis complete!", type = "message")
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
  
  observeEvent(input$toggleSignatureSidebar, {
    toggle("signatureSidebar", anim = TRUE, animType = "slide")
  })
  
  observeEvent(input$toggleCorrelationSidebar, {
    toggle("correlationSidebar", anim = TRUE, animType = "slide")
  })
  
  observeEvent(input$toggleSurvivalSidebar, {
    toggle("survivalSidebar", anim = TRUE, animType = "slide")
  })

  output$signatureCountText <- renderText({
    req(input$pathway)
    total <- length(unique(merged_sing_df$Pathway))
    paste("Selected", length(input$pathway), "of", total, "Signatures")
  })
  },
   error = function(e) {
    print(paste("SERVER ERROR:", e$message))
  })
}