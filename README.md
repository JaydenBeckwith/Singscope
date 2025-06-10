# Singscope

<p align="center">
  <img src="https://github.com/user-attachments/assets/4b44ca28-5031-4d1b-b5a2-ef057ebf8b44" alt="Singscope logo" width="400"/>
</p>

An interactive Shiny application for analysing bulk RNA signatures and supporting clinical stratification of melanoma patient cohorts.
---

## Features

<p align="center">
  <img src="https://github.com/user-attachments/assets/9fb8c652-f99c-43f5-b12f-b87fcb5c62a1" alt="toolexplanation" width="800"/>
</p>

* Upload and visualise signature enrichment scores (Singscores)
* Track treatment response and recurrence distributions
* Correlation heatmaps of selected signatures
* Track pathway correlations over time
* Survival analysis with KM plots
---
## Installation (Docker)

1. **Clone this repo**

   ```bash
   git clone https://github.com/JaydenBeckwith/Singscope.git
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
   cd /srv/shiny-server/R
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
