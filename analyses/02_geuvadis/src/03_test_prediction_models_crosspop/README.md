# Train/test prediction models: GEUVADIS constituent pops

Use elastic net to train and test predictive models with GEUVADIS data.
Unlike the [previous analysis step](../02_test_prediction_models_continentalpop/README.md), this analysis partitions samples into five constituent 1000 Genomes populations:
* CEU, the Central Europeans from Utah
* GBR, British from Great Britain
* TSI, Tuscans from Italy
* FIN, Finnish from Finland
* YRI, Yoruba from Nigeria

All samples are unrelated; the children in the GEUVADIS dataset are excluded.
Each population is subsampled to the size of the smallest population (YRI, _n_ = 89).
Prediction models are trained in one population and tested in all five populations.

## Prerequisites
* `Rscript`
* `qsub`

The following R packages are used:
* `glmnet`
* `dplyr`
* `data.table`
* `broom`
* `optparse`
* `dunn.test`
* `ggplot2`

Install there from [CRAN](https://cran.r-project.org) using any standard installation procedure. One way is to type
```R
install.packages(c("glmnet", "data.table", "dplyr", "optparse", "broom", "dunn.test", "ggplot2"))
```
from within the R console.

## Running 

From the command line in an SGE environment, execute
```bash
./test_prediction_models_crosspop.sh
```
This will schedule jobs with the correct ordering and hold patterns.
Jobs for each train-test case are broken into prediction, collection, and postprocessing.
Since each training population is independent, all five training populations can be run simultaneously.

Once all results are available, they can be analyzed by running
```bash
./run_geuvadis_compile_prediction_results_onepop.sh
```
which will compile results for each training population individually and then concatenate them together.

## Notes

Similar to the previous analysis step, this analysis requires substantial disc space (> 100Gb).
Delete log files as needed once correct execution is confirmed.
