


run_survival_analysis <- function(data, survival_type = c("OS", "RFS", "EFS", "MSS"),
                                  group_col, time_unit = c("days", "months", "years")) {
  
  
  survival_type <- match.arg(survival_type)
  time_unit <- match.arg(time_unit)

  data <- data %>%
    dplyr::mutate(across(
      c(
        date_neoadj_start, date_surgery, recurrencefree_date, eventfree_date, date_fu,
        recurrence_status, event_status, death_cause
      ),
      ~ dplyr::na_if(as.character(.), ".")
    ))


  ### make sure we are only counting each patient once - since we can have multiple samples per ptx
  data <- data %>% dplyr::distinct(patient_id, .keep_all = TRUE)

  # Convert Excel-style numeric dates or string dates
  date_cols <- c(
    "date_fu", "date_neoadj_start", "date_surgery",
    "recurrencefree_date", "eventfree_date"
  )


  ### R uses american dates - correct for this 
  for (col in date_cols) {
    if (is.numeric(data[[col]])) {
      data[[col]] <- as.Date(data[[col]], format = "%d-%m-%Y")
    } else {
      data[[col]] <- as.Date(as.character(data[[col]]), format = "%d/%m/%Y")
    }
  }

  # Ensure key columns are correctly formatted
  data <- data %>%
  dplyr::mutate(
    recurrence_status = as.numeric(recurrence_status),
    event_status = as.numeric(event_status),
    death = as.numeric(death),
    death_cause = as.numeric(death_cause),
    surgery = as.numeric(surgery)
  )

# Define eligibility per survival type as per surv rules 
data <- data %>%
  dplyr::mutate(
    eligible_OS  = !is.na(date_neoadj_start) & !is.na(date_fu),
    eligible_EFS = !is.na(date_neoadj_start) & !is.na(eventfree_date),
    eligible_RFS = surgery == 1 & !is.na(date_surgery) & !is.na(recurrencefree_date) & !is.na(recurrence_status),
    eligible_MSS = !is.na(date_neoadj_start) & !is.na(date_fu) & !is.na(death_cause)
  ) 


# Define raw time and event columns using proper logic
data <- data %>%
  dplyr::mutate(
    raw_time = dplyr::case_when(
      survival_type == "OS"  & eligible_OS  ~ as.numeric(date_fu - date_neoadj_start),
      survival_type == "EFS" & eligible_EFS ~ as.numeric(eventfree_date - date_neoadj_start),
      survival_type == "RFS" & eligible_RFS ~ as.numeric(recurrencefree_date - date_surgery),
      survival_type == "MSS" & eligible_MSS ~ as.numeric(date_fu - date_neoadj_start),
      TRUE ~ NA_real_
    ),
    event = dplyr::case_when(
      survival_type == "OS"  ~ as.numeric(death == 1),
      survival_type == "EFS" ~ as.numeric(event_status == 1),
      survival_type == "RFS" ~ as.numeric(recurrence_status == 1),
      survival_type == "MSS" ~ as.numeric(death_cause == 1),
      TRUE ~ NA_real_
    ),
    time_unit_factor = dplyr::case_when(
      time_unit == "days"   ~ 1,
      time_unit == "months" ~ 1 / 30.44,
      time_unit == "years"  ~ 1 / 365.25
    ),
    time = raw_time * time_unit_factor
  )

print(data)

# Filter incomplete or malformed rows
#data <- data %>%
  #filter(!is.na(time), !is.na(event), time >= 0, !is.na(.data[[group_col]]))

# Stop if nothing usable
if (nrow(data) == 0) {
  stop("No valid observations for survival analysis after filtering.")
}

  # Fit survival model
  # Create formula
  surv_formula <- as.formula(paste("survival::Surv(time, event) ~", group_col))

  # Fit survival model
  surv_fit <- survival::survfit(surv_formula, data = data)
  surv_fit$call$formula <- surv_formula
  surv_diff <- survival::survdiff(surv_formula, data = data)

  if (!is.null(surv_fit$strata)) {
  names(surv_fit$strata) <- gsub(".*=", "", names(surv_fit$strata))
  surv_fit$strata <- setNames(surv_fit$strata, names(surv_fit$strata))
}


  pval <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  pval_txt <- paste0("Log-rank p = ", signif(pval, 3))

  surv_data <- list(
    fit = surv_fit,
    data = data,
    pval_txt = pval_txt,
    group_col = group_col,
    time_unit = time_unit
  )

  return(surv_data)
}

# === Render Plot Wrapper ===
plotly_survival <- function(surv_fit, data, survival_type, group_col, pval_txt, time_unit = "months") {

  group_col <- as.character(group_col)

  data[[group_col]] <- as.factor(data[[group_col]])
  levels(data[[group_col]]) <- gsub(".*=", "", levels(data[[group_col]]))

  #set max on plot so it doesnt hang over 

  xmax <- quantile(data$time, 0.95, na.rm = TRUE)

  ggs <- survminer::ggsurvplot(
    fit = surv_fit,
    data = data,
    risk.table = TRUE,
    ncensor.plot = TRUE,
    pval = FALSE,
    legend.title = group_col,
    xlab = time_unit,
    ylab = paste("Survival Probability -", survival_type),
    censor.shape = "|",
    censor.size = 3,
    risk.table.y.text.col = TRUE,
    palette = "Set1",
    conf.int.style = "step",
    surv.median.line = "hv",
    ggtheme = theme_minimal(),
    xlim = c(0, xmax)
  )

  # Main plot adjustments
  ggs$plot <- ggs$plot + labs(title = pval_txt) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(size = 18, face = "bold")
    )

  plotly_plot <- plotly::ggplotly(ggs$plot)

  # Format risk table & censor plot
  ggs$table <- ggs$table +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )

  ggs$ncensor.plot <- ggs$ncensor.plot + theme(
      legend.position = "bottom",
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )

  return(list(
    plot = plotly_plot,
    gplot = ggs$plot,     # the original ggplot object
    risk_table = ggs$table,
    censor_plot = ggs$ncensor.plot
  ))
}