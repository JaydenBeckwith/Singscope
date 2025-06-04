library(survival)
library(survminer)


run_survival_analysis <-  function(data, survival_type = c("OS", "RFS"), group_col = "Response", title = NULL) {
  survival_type <- match.arg(survival_type)

  # Determine column names
  time_col <- if (survival_type == "OS") "os_time" else "rfs_time"
  event_col <- if (survival_type == "OS") "os_status" else "rfs_status"
  if (is.null(title)) {
    title <- if (survival_type == "OS") "Overall Survival" else "Recurrence-Free Survival"
  }

  # Validate presence of required columns
  required <- c(time_col, event_col, group_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

  surv_obj <- Surv(time = data[[time_col]], event = data[[event_col]])
  fit <- survfit(surv_obj ~ data[[group_col]])

  surv_diff <- survdiff(surv_obj ~ data[[group_col]])
  pval <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  pval_txt <- paste0("Log-rank p = ", signif(pval, 3))

  # Plot with ggsurvplot
  ggs <- ggsurvplot(
    fit,
    data = data,
    risk.table = TRUE,
    pval = FALSE,
    legend.title = group_col,
    ggtheme = theme_minimal(),
    title = title
  )

  # Convert to plotly for Shiny compatibility
  ggplotly(ggs$plot) %>%
    layout(title = list(text = paste0(title, "<br><sup>", pval_txt, "</sup>")))
}