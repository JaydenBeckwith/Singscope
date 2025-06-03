# Singscope

<p align="center">
  <img src="https://github.com/user-attachments/assets/c0c183eb-d235-4fc2-a005-cc6c4fe31a53" alt="Singscope logo" width="300"/>
</p>

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
docker run --rm --workdir /srv/shiny-server/tests singscope Rscript run_tests.R
```

---
## Singscore Citation 

- Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ (2018).  
  *Single sample scoring of molecular phenotypes.*  
  **BMC Bioinformatics**, 19(1), 404.  
  [https://doi.org/10.1186/s12859-018-2435-4](https://doi.org/10.1186/s12859-018-2435-4)

- Bhuva DD, Cursons J, Davis MJ (2020).  
  *Stable gene expression for normalisation and single-sample scoring.*  
  **Nucleic Acids Research**, 48(19), e113.  
  [https://doi.org/10.1093/nar/gkaa802](https://doi.org/10.1093/nar/gkaa802)
---