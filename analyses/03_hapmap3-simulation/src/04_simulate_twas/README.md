# Simulate transcriptome-wide association analysis 
This analysis uses output from the [previous analysis](../03_compile_results/README.md).

## Analysis steps
There are two steps, which must occur in a particular order:
1. Simulate the TWAS using simulated and predicted gene expression. 
2. Plot results 

Step (1) schedules parallel jobs similar to the [gene expression simulations](../02_test_prediction_models/README.md).
Because this is a serial loop, it can take a while to complete.
To avoid timeouts, step (1) can be run in the background.
`screen` accomplishes this neatly. So does `nohup`:
```bash
nohup Rscript compile_all_results.R \
    --results-directory ${output_data_dir} \
    --output-directory ${resultsdir} \
    --output-filename ${results_file} \
    > nohup.compile.results.out 2> nohup.compile.results.err &
```

Step (2) saves plots directly to file, but outputs results from statistical analysis to the console.
As in the example above, the results can be redirected and saved to file:
```bash
Rscript plot_results.R \
     --results-file ${results_filepath} \
     --output-directory ${plotdir} \
     --plot-filetype ${plot_filetype} \
    > plot.results.out 2> plot.results.err
```

## Running
From a Bourne shell, call the script directly:
```bash
./simulate_twas.sh
```
and then compile results and plots with
```bash
Rscript plot_twas_results.R
```
