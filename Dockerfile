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
COPY packages.txt /tmp/packages.txt

RUN Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org')"

# Install specific versions
COPY packages.txt /tmp/packages.txt
RUN Rscript -e "\
  pkgs <- readLines('/tmp/packages.txt'); \
  for (pkg in pkgs) { \
    parts <- strsplit(pkg, '==')[[1]]; \
    name <- parts[1]; \
    version <- parts[2]; \
    remotes::install_version(name, version = version, repos = 'https://cloud.r-project.org') \
  }"

RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" && \
    Rscript -e "BiocManager::install(version = '3.18')" && \
    Rscript -e "BiocManager::install(c('biomaRt', 'SummarizedExperiment', 'singscore'))"

# Copy your Shiny app to the image
COPY . /srv/shiny-server/
RUN chown -R shiny:shiny /srv/shiny-server

# Expose the Shiny port
EXPOSE 3838

# Run the Shiny server
CMD ["shiny-server"]
