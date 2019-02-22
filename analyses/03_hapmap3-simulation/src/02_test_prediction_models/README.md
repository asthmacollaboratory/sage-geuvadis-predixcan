# Test crosspopulation prediction

This folder contains code to construct predictive models and test them in various populations.
Many model scenarios are tested. The relevant parameters are

| Variable | Parameter | Comments |
| --- | --- | --- |
| `model_sizes` | number(s) of causal eQTLs in a model | bounded by maximum gene length, but for analysis purposes varies from 1 to 20 |
| `seeds` | random seeds | must have at least 1 to make random sampling reproducible |
| `props` | shared eQTL proportions | must be between 0 and 1 |
| `same_eqtls` | do eQTL positions match? | either `TRUE` or `FALSE`. Setting to `TRUE` overrides information about shared eQTL proportions |
| `same_effects` | do eQTL effects match? | either `TRUE` or `FALSE`, though cases with `FALSE` are not considered here. |


## Analysis steps
1. Pull top genes from list of gene positions downloaded in [analysis 1](../01_simulate_populations/README.md).
2. Write each testing scenario as a job to be read from a master file. 
3. Schedule and execute jobs in parallel. 

Step (2) uses elastic net regression from the `glmnet` package to estimate predictive models in a nested cross-validation scheme.
The function `cv.glmnet` runs the internal cross-validation loop, while the external folds are partitioned manually.
Four variables control the cross-validation behavior:

| Variable | Parameter | Comments |
| --- | --- | --- |
| `nfolds_external` | number of external cross-validation folds | default: 5 |
| `nolds_internal` | number of internal cross-validation folds | controlled by `glmnet`, default: 10 |
| `nfolds_parallel` | number of internal folds run in parallel | also controlled by `glmnet`, default: 1 | 
| `nthreads` | number of concurrent threads | one thread per job, default: 24 |

It is _strongly recommended_ to keep `nfolds_parallel=1` since the internal folds are computationally light.
Instead, maximize throughput by parallelizing across jobs (testing scenarios) with `nthreads`.

Step (3) relies heavily on [`nohup`](https://en.wikipedia.org/wiki/Nohup) and schedules jobs to run in the background. 
In principle, this analysis can use all available virtual CPU cores, but in practice `libgomp` may limit the number of jobs that can run in parallel.
The manuscript used `nthreads=24`.
Your mileage may vary.

Note that step (3) is easily amenable to running with `qsub` on an SGE cluster if the user has access to one.
The commands written to `simulation_joblist.sh` can be prefaced with `qsub` and its requisite options when writing `simulation_joblist.sh`.
Afterwards, call `simulation_joblist.sh` directly to place all jobs in a queue.

## Running
From a Bourne shell, call the script directly:
```bash
./test_prediction_models.sh
```


## References
* Gamazon, E.R., Wheeler, H.E., Shah, K.P., Mozaffari, S.V., Aquino-Michaels, K., Carroll, R.J., Eyler, A.E., Denny, J.C., Nicolae, D.L., Cox, N.J., _et al._ (2015) A gene-based association method for mapping traits using reference transcriptome data. _Nature Genetics_ 47, 1091–1098.
* Friedman J, Hastie T, and Tibshirani R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. _Journal of Statistical Software_, **33**(1), 1-22. ([link](http://www.jstatsoft.org/v33/i01/.))
* Zou, H. and Hastie, T. (2005) Regularization and variable selection via the elastic net. _Journal of the Royal Statistical Society Series B_ **67**(2), 301-320.
