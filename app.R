# Load necessary libraries
if (!require(shiny)) install.packages("shiny")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(DT)) install.packages("DT")

options(shiny.maxRequestSize = 400*1024^2)  # 400 MB
source("preprocess_data.R")
source("ui.R")
source("server.R")

# === Run the Application ===
shinyApp(ui = ui, server = server)