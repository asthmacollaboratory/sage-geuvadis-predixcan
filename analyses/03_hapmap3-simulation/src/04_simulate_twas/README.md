# Simulate transcriptome-wide association analysis 
This analysis uses output from the [previous analysis](../03_compile_results/README.md).

## Analysis steps
There are two steps, which must occur in a particular order:
1. Simulate the TWAS using simulated and predicted gene expression. 
2. Plot results 

Step (1) schedules parallel jobs similar to the [gene expression simulations](../02_test_prediction_models/README.md).
The plotting in step (2) must be run manually after step (1) finishes.

## Running
From a Bourne shell, call the script directly:
```bash
./simulate_twas.sh
```
and then compile results and plots with
```bash
Rscript plot_twas_results.R
```
