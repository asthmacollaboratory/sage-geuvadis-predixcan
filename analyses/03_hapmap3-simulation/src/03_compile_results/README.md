# Compile prediction results
This analysis compiles the output from the [second analysis](../02_test_prediction_models/README.md).

## Analysis steps
There are only two steps here:
1. Compile the data into a single `data.table` and write to file
2. Run statistical analyses and create plots

Step (1) loops through the output from the second analysis.
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
./compile_results.sh
```

## References
* Dinno A. (2017) `dunn.test`: Dunn's Test of Multiple Comparisons Using Rank Sums. R package version 1.3.5. ([link](https://CRAN.R-project.org/package=dunn.test))
* Wickham H. (2016) `ggplot2`: Elegant Graphics for Data Analysis. Springer-Verlag New York.
