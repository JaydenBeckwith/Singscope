# Start from Rocker Shiny base image with R 4.3
FROM rocker/shiny:4.3.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN R -e "install.packages(c( \
    'shiny', 'shinythemes', 'shinyjs', 'shinyWidgets', 'plotly', 'DT', \
    'shinycssloaders', 'dplyr', 'tidyr', 'tibble', 'readr', 'RColorBrewer', \
    'tidyverse' \
), repos = 'https://cloud.r-project.org')"

# Install Bioconductor and packages
RUN R -e \"if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\" \\\n    && R -e \"BiocManager::install(version = '3.18')\" \\\n    && R -e \"BiocManager::install(c('biomaRt', 'SummarizedExperiment'))\"

# Copy your Shiny app to the image
COPY . /srv/shiny-server/
RUN chown -R shiny:shiny /srv/shiny-server

# Expose the Shiny port
EXPOSE 3838

# Run the Shiny server
CMD [\"/usr/bin/shiny-server\"]