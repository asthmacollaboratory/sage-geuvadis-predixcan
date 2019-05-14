#!/usr/bin/env bash
#
# DO NOT RUN THIS SCRIPT DIRECTLY
# The code below is offered for demonstration purposes only

# run first step
./01_simulate_populations/./simulate_populations.sh

# run second step
./02_test_prediction_models/./test_prediction_models.sh

# run third step
./03_compile_results/./compile_results.sh

# run fourth step
./04_simulate_twas/./simulate_twas.sh
