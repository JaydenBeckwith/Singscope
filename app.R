# Load necessary libraries
if (!require(shiny)) install.packages("shiny")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")
if (!require(DT)) install.packages("DT")

library(dplyr)
library(shiny)
library(tidyr)
library(DT)
library(RColorBrewer)

options(shiny.maxRequestSize = 400*1024^2)  # 400 MB
source("preprocess_data.R")
source("ui.R")
source("server.R")

### run 
shinyApp(ui = ui, server = server)
