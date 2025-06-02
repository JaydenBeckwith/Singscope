# Start from Rocker Shiny base image with R 4.3
FROM rocker/shiny:4.3.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
  libglpk40 \
  libglpk-dev \
  libxml2-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libxt-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  zlib1g-dev \
  build-essential

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
    Rscript -e "BiocManager::install(c('biomaRt', 'SummarizedExperiment', 'singscore'))" && \
    Rscript -e "BiocManager::install('clusterProfiler', ask = FALSE, force = TRUE)"

# Copy your Shiny app to the image
COPY . /srv/shiny-server/
RUN chown -R shiny:shiny /srv/shiny-server

# Expose the Shiny port
EXPOSE 3838

# Run the Shiny server
CMD ["shiny-server"]
