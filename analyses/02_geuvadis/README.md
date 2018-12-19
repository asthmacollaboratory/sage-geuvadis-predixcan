# Crosspopulation prediction analysis using GEUVADIS genotype and gene expression data 
Pipeline to train and analyze [PrediXcan](https://github.com/hakyim/PrediXcan) prediction weights using data from the Genetic European Variation in Health and Disease ([GEUVADIS](www.geuvadis.org/web/geuvadis)) project.

Genotypes and normalized gene expression measures for 373 EUR samples and 89 YRI samples from GEUVADIS were downloaded from the [EBI website](https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/analysis_results/).

Per protocol from [Gamazon _et al._ (2015)](https://www.ncbi.nlm.nih.gov/pubmed/26258848), new prediction weights for each population were computed using [elastic net regularized](https://en.wikipedia.org/wiki/Elastic_net_regularization) linear regression with the [`glmnet`](https://cran.r-project.org/web/packages/glmnet/index.html) package in R.

## Steps
1. Download and parse GEUVADIS data ([docs](./src/01_download_parse_data/README.md))
