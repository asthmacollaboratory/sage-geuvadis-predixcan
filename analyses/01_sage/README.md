# Analysis of Predixcan predictions in SAGE 
This repository documents code used in the comparison of PrediXcan predictions versus RNA-Seq measurements using pilot data from the Study of African Americans, Asthma, Genes, and Environment (SAGE).

## Data

Genotype data are available on dbGaP under accession number phs000921 (see [here](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000921.v3.p1)).
Expression data are available upon request by contacting the [UCSF Asthma Collaboratory](https://pharm.ucsf.edu/burchard/contact).

PrediXcan prediction weights were pulled from [PredictDB](http://predictdb.org/).

## Analysis steps
1. Format genotype and gene expression data ([docs](./src/01_parse_data/README.md))
2. Run PrediXcan with genotype data to predict gene expression ([docs](./src/02_run_predixcan/README.md))
3. Compile prediction results and compare with measurements ([docs](./src/03_analyze_results/README.md))

## References
* Borrell LN, _et al_. Childhood obesity and asthma control in the GALA II and SAGE II studies. _AJRCCM_ 2013 Apr 01; 187(7):697-702. ([link](https://www.ncbi.nlm.nih.gov/pubmed/23392439))
* Mak ACY, White MJ, Eckalbar W, _et al_. Whole Genome Sequencing Study on Bronchodilator Drug Response in Ethnically Diverse Children with Asthma. _AJRCCM_ 2018; 197(12):1552--1564. ([link](https://www.ncbi.nlm.nih.gov/pubmed/29509491)
* Gamazon E, _et al_. A gene-based association method for mapping traits using reference transcriptome data. _Nature Genetics_ 2015; 47:1091--1098. ([link](https://www.nature.com/articles/ng.3367))
