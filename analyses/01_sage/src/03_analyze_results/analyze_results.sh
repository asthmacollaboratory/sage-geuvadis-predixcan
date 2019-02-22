#!/usr/bin/env bash
# ==========================================================================================
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys
# ==========================================================================================

# ==========================================================================================
# script variables 
# ==========================================================================================

# directories
thisdir="$(dirname $(readlink -f $0))"

analysisdir="${thisdir}/../../analysis"
datafiles_dir="${thisdir}/../../datafiles"

datadir="${analysisdir}/data"

rnaseqdir="${datadir}/rnaseq"

resultsdir="${analysisdir}/results"
predixcanoutdir="${resultsdir}/predixcan"

# change plot type depending on publication avenue
# default is PNG
plot_type="png"

$RSCRIPT $R_analyze_predictions \
    --predixcan-dir ${predixcanoutdir} \
    --data-dir ${rnaseqdir} \
    --output-dir ${resultsdir} \
    --plot-type ${plot_type}
