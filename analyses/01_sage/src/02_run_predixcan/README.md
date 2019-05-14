# Run PrediXcan on SAGE data

This analysis step runs PrediXcan on SAGE data and postprocesses the output for analysis.

# Running
In a BASH console in an SGE environment, type
```bash
./run_predixcan.sh
```

## Notes

It is wise to run PrediXcan using _imputed_ genotypes.
For this study, genotypes were imputed on the [Michigan Imputation Server](https://imputationserver.sph.umich.edu) with the following settings:
* phasing: EAGLE
* reference panel: 1000 Genomes phase 3 version 5
