# Train/test prediction models: GEUVADIS continental populations

Use elastic net to train and test predictive models with GEUVADIS data.
The population breakdown is as follows:

| Population | Description | Training pop? | Testing pop? |
| --- | --- | --- | --- |
| EUR373 | 373 individuals from four 1000 Genomes populations: CEU + GBR + TSI + FIN | yes | yes |
| EUR278 | 278 individuals from three 1000 Genomes populations: CEU + GBR + TSI | yes | yes |
| FIN | 95 individuals from the 1000 Genomes Finnish population (FIN) | no | yes |
| AFR | 89 individuals from the 1000 Genomes Yoruba population (YRI) | yes | yes | 

All samples are unrelated; the children in the GEUVADIS dataset are excluded.

## Prerequisites
* `Rscript`
* `qsub`

The following R packages are used:
* `glmnet`
* `dplyr`
* `data.table`
* `broom`
* `optparse`

Install there from [CRAN](https://cran.r-project.org) using any standard installation procedure. One way is to type
```R
install.packages(c("glmnet", "data.table", "dplyr", "optparse", "broom"))
```
from within the R console.

## Running 

From the command line in an SGE environment, execute
```bash
./test_prediction_models_continentalpop.sh
```
This will schedule jobs with the correct ordering and hold patterns. Jobs for each train-test case are broken into prediction, collection, and postprocessing.

## Notes

As coded, the case of EUR278 to FIN cannot run in parallel with EUR278 to AFR, since it creates race conditions on the output files. Schedule it once all other jobs are completed. 
Beware that between output files and log files from job execution, this analysis requires substantial storage space (> 100Gb).
Jobs will fail (possibly quietly) once disc capacity is reached.
The log files can be safely deleted after checking for correct execution.
