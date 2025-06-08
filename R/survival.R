


run_survival_analysis <- function(data, survival_type = c("OS", "RFS", "EFS", "MSS"), group_col) {
  survival_type <- match.arg(survival_type)

  data <- data %>%
  dplyr::mutate(across(
    c(
      date_neoadj_start, date_surgery, recurrencefree_date, eventfree_date, date_fu,
      recurrence_status, event_status, death_cause
    ),
    ~ dplyr::na_if(as.character(.), ".")
  ))

   # Convert Excel-style date columns
  date_cols <- c(
    "date_fu", "date_neoadj_start", "date_surgery",
    "recurrencefree_date", "eventfree_date"
  )
  
  for (col in date_cols) {
  if (is.numeric(data[[col]])) {
    data[[col]] <- as.Date(data[[col]], origin = "1899-12-30")
  } else {
    # Try parsing string dates like "2020-01-01"
    data[[col]] <- as.Date(data[[col]])
  }
}

  # Ensure key columns are in correct format
  data <- data %>%
    dplyr::mutate(
      recurrence_status = as.numeric(recurrence_status),
      event_status = as.numeric(event_status)
    )

  # Compute time and event variables
  data <- data %>%
  dplyr::mutate(
    time = dplyr::case_when(
      survival_type == "OS"  ~ as.numeric(as.Date(date_fu) - as.Date(date_neoadj_start)),
      survival_type == "RFS" ~ as.numeric(as.Date(recurrencefree_date) - as.Date(date_surgery)),
      survival_type == "EFS" ~ as.numeric(as.Date(eventfree_date) - as.Date(date_neoadj_start)),
      survival_type == "MSS" ~ as.numeric(as.Date(date_fu) - as.Date(date_neoadj_start))
    ) / 30.44,  # Convert to months
    event = dplyr::case_when(
      survival_type == "OS"  ~ as.numeric(death == 1),
      survival_type == "RFS" ~ as.numeric(recurrence_status == 1),
      survival_type == "EFS" ~ as.numeric(event_status == 1),
      survival_type == "MSS" ~ as.numeric(death_cause == 1)
    )
  )

  # Apply RFS-specific eligibility rules
  if (survival_type == "RFS") {
    data <- data %>%
      filter(surgery == 1, !is.na(date_surgery), !is.na(recurrencefree_date), !is.na(recurrence_status))
  }

  # Filter incomplete/malformed rows
  data <- data %>%
    filter(!is.na(time), !is.na(event), time >= 0, !is.na(.data[[group_col]]))

  # Abort early if no data left
  if (nrow(data) == 0) {
    stop("No valid observations for survival analysis after filtering.")
  }

  # Create formula
  surv_formula <- as.formula(paste("survival::Surv(time, event) ~", group_col))

  # Fit survival model
  surv_fit <- survival::survfit(surv_formula, data = data)

  surv_fit$call$formula <- surv_formula

  surv_diff <- survival::survdiff(surv_formula, data = data)

  # Compute log-rank p-value
  pval <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  pval_txt <- paste0("Log-rank p = ", signif(pval, 3))
  

  surv_data = list(
    fit = surv_fit,
    data = data,
    pval_txt = pval_txt,
    group_col = group_col
  )

  return(surv_data)

}



# === Render Plot Wrapper ===
plotly_survival <- function(surv_fit, data, survival_type, group_col, pval_txt) {
  group_col <- as.character(group_col)

  ggs <- survminer::ggsurvplot(
    fit = surv_fit,
    data = data,
    risk.table = TRUE,
    pval = pval_txt,
    legend.title = group_col,
    xlab = "Months",
    ylab = paste("Survival Probability -", survival_type),
    censor.shape = "|",
    censor.size = 3,
    palette = "Set1",
    ggtheme = theme_minimal()
  )

  return(list(
    plot = ggs$plot,
    risk_table = ggs$table
  ))
}
