# Download and parse GEUVADIS data

This analysis step compiles GEUVADIS genotype and gene expression data for analysis

## Prerequisites
* an internet connection (scripts here download from the EBI server)
* `bcftools`
* `plink`
* `bgzip`
* `Rscript`
* `wget`
* `gzip`

R packages used here include:
* `data.table`
* `dplyr`
* `optparse`

Install there from [CRAN](https://cran.r-project.org) using any standard installation procedure. One way is to type
```R
install.packages(c("data.table", "dplyr", "optparse"))
```
from within the R console.

## Running 

Genotype and expression data each have their own script and can be run in parallel. Execute normally, or use `nohup` or `screen` to run them in the background:
```bash
./download_parse_genotype_data.sh
./download_parse_expression_data.sh
```

## Notes

The total space occupied by GEUVADIS VCFs is reasonably large. The analysis should be run with at least 120 gigabytes of available disc space.
Note that a PLINK binary genotype file is parsed from the VCFs, so the VCFs can be deleted once the PLINK file is available. 
Deleting VCFs will recover approximately 100-110Gb of disc space. 
