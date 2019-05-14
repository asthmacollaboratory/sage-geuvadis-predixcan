# Crosspopulation prediction analysis using HapMap3 data
Analysis of crosspopulation prediction capacity using a simulation of two ancestral populations and an admixed one from HapMap3 phase 3 data.

## Prerequisites
* Python 2.7 (note: scripts here are incompatible with Python 3) 
* [R](https://www.r-project.org/) version 3.4.3 or higher
* a Linux computing environment

All other requisite binaries and external data are downloaded and configured automatically.

This analysis makes extensive use of `Rscript`. 
It assumes a correct installation and configuration of the default R environment. 

The R packages used here are
* `tidyverse`, particularly `ggplot2` and `dplyr` 
* `knitr`
* `data.table`
* `optparse`
* `doParallel`
* `glmnet`
* `assertthat`

All of these packages are registered on [CRAN](https://cran.r-project.org/).
Install these packages with any [standard approach](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).
From within R, a one-shot approach to installing packages is to type 
```R
install.packages(c(tidyverse, knitr, data.table, optparse, doParallel, glmnet, assertthat))
```
assuming compilers, library paths, CRAN repository, and write permissions are all configured correctly. 

Additionally, the following Bioconductor packages are required:
* [`biomaRt`](https://www.bioconductor.org/packages/release/bioc/html/biomaRt.html)
* [`annotate`](https://www.bioconductor.org/packages/release/bioc/html/annotate.html)

Install these packages using the [Bioconductor protocol](https://www.bioconductor.org/install/).

## Running
The analysis is divded into three parts, each with its own documentation:
1. Simulating haplotypes ([docs](./src/01_simulate_populations/README.md))
2. Estimating and evaluating prediction models ([docs](./src/02_test_prediction_models/README.md))
3. Compiling and plotting results ([docs](./src/03_compile_results/README.md))
4. Simulating a TWAS ([docs](./src/04_simulate_twas/README.md))

The analysis will generate a directory `./src/analysis` with all output files.

Each analysis step has its own BASH script. Run the scripts in sequential order. A demonstration script is provided for this purpose:
```bash
./src/./run_analysis.sh
```
Note that step (2) uses a parallel computing framework that sends jobs into the shell background. The user must wait for all background jobs to complete before proceeding to step (3).


## Notes
Simulations from step (1) are not reproducible. Simulated data used in the full analysis can be found [here](https://ucsf.box.com/s/fgjizb0uiopr6sw76tf0u5beg7js2l36).

The parameters for testing different prediction scenarios in step (2) are set for demonstration purposes only.

The full analysis tests
- 4 model sizes
- 100 random seeds
- 98 genes
- 11 shared eQTL proportions
for a total of 431,200 jobs. Consequently, the full analysis generates a **HUGE** amount of data.
**DO NOT RUN** a full analysis without at least 1.5TB of disc space!
