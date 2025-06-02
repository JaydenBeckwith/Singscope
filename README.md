# Singscope

An interactive Shiny application for exploring pathway-level singscore enrichment data, correlation heatmaps, and clinical stratification for melanoma cohorts.
---

## Features

* Upload and visualise signature enrichment scores
* Compare mutation and response distributions
* Correlation heatmaps of selected pathways
* Track pathway correlations over time
* Survival analysis with KM plots
---
## Installation (Docker)

1. **Clone this repo**

   ```bash
   git clone https://github.com/yourusername/Singscope.git
   cd Singscope
   ```

2. **Build the Docker image**

   ```bash
   docker build -t singscope .
   ```

3. **Run the app locally**

   ```bash
   docker run -p 3838:3838 singscope
   ```

4. **Open in your browser**

   ```
   http://localhost:3838
   ```

## Run interactive Docker for logging 

```bash
   docker run -it --rm -p 3838:3838 singscope bash
   cd /srv/shiny-server
   R -e "shiny::runApp('.', port = 3838, host = '0.0.0.0')"
```

---

## Run Tests

You can run `testthat` tests directly:

```bash
docker run --rm singscope Rscript -e "testthat::test_dir('tests/testthat')"
```

---