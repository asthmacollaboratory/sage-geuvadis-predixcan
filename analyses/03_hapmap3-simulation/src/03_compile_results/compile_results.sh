#!/usr/bin/env bash
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
# ==========================================================================================


# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -e  ## script will exit on error
set -u  ## script will exit if it sees an uninitialized variable


# ==========================================================================================
# binaries, directories, filepaths 
# ==========================================================================================
RSCRIPT=$(whereis Rscript | awk '{print $2}') ## will default to system Rscript, if installed 

thisdir="$(dirname $(readlink -f $0))"
analysisdir="${thisdir}/../../analysis"
outputdir="${analysisdir}/prediction_output"
output_text_dir="${outputdir}/prediction_output_text"
output_data_dir="${outputdir}/prediction_output_data"
plotdir="${analysisdir}/plots"
resultsdir="${analysisdir}/results"


# ==========================================================================================
# external scripts 
# ==========================================================================================
R_compile_results="${thisdir}/compile_all_results.R"
R_plot_results="${thisdir}/plot_results.R"

# ==========================================================================================
# script variables
# ==========================================================================================
results_file="1kg.all.results.txt"
results_filepath="${resultsdir}/${results_file}"
plot_filetype="pdf"

# ==========================================================================================
# executable code 
# ==========================================================================================

# make directories if they don't exist
mkdir -p ${plotdir}
mkdir -p ${resultsdir}

# compile results and save to $resultsdir
$RSCRIPT $R_compile_results \
    --results-directory ${output_data_dir} \
    --output-directory ${resultsdir} \
    --output-filename ${results_file}

# plot results and save to $plotdir
$RSCRIPT $R_plot_results \
    --results-file ${results_filepath} \
    --output-directory ${plotdir} \
    --plot-filetype ${plot_filetype}
