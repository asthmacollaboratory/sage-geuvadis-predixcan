#!/usr/bin/env bash
# ==========================================================================================
# copyright Asthma Collboratory (2018)
# coded by Kevin L. Keys
#
# This script computes prediction weights from GEUVADIS transcriptome data.
# The prediction models are estimated in accordance with PrediXcan protocols;
# for reference, see:
# -- Gamazon et al. (2015) Nature Genetics
# -- https://github.com/hakyimlab/PrediXcan/tree/master/Paper-Scripts
#
# The GEUVADIS populations used here are:
# 1. EUR373, which includes all 373 unrelated Europeans from GEUVADIS
# 2. EUR278, which includes all 278 non-Finnish Europeans
# 3. FIN, which includes all 95 Finns from GEUVADIS
# 4. YRI, which includes all 89 Yoruba from GEUVADIS
#
# This script schedules four train-test scenarios:
# -- EUR373 to YRI
# -- YRI to EUR373
# -- EUR278 to YRI
# -- EUR278 to FIN
#
# Call:
#
# ./test_prediction_models_continentalpop.sh
# ==========================================================================================


# ==========================================================================================
# BASH script settings
# ==========================================================================================

set -o errexit # set -e, script will exit on error
set -o nounset # set -u, script will exit if it sees an uninitialized variable
#set -o xtrace  # set -x, script will track which command is currently running

# ==========================================================================================
# script variables
# ==========================================================================================

thisdir="$(dirname $(readlink -f $0))"
analysisdir="${thisdir}../analysis"
resultsdir="${analysisdir}/results/glmnet/elasticnet"

BASH_schedule_eur373_to_yri="${thisdir}/schedule_eur373_to_yri.sh"
BASH_schedule_yri_to_eur373="${thisdir}/schedule_yri_to_eur373.sh"
BASH_schedule_eur278_to_fin="${thisdir}/schedule_eur278_to_fin.sh"
BASH_schedule_eur278_to_yri="${thisdir}/schedule_eur278_to_yri.sh"
BASH_define_variables="${thisdir}/../common/geuvadis_variables.sh"
BASH_compile_results_continentalpop="${thisdir}/qsub_geuvadis_compile_continentalpop_results.sh"

R_compile_results_continentalpop="${thisdir}/geuvadis_compile_largepop_results.R"


# ==========================================================================================
# source variables
# ==========================================================================================

# must first grab all of the variables for GEUVADIS analysis
source ${BASH_define_variables}


# ==========================================================================================
# start job scheduling script
# ==========================================================================================

# -------------------- #
# EUR373 to YRI

source ${BASH_schedule_eur373_to_yri}

# -------------------- #
# YRI to EUR373

source ${BASH_schedule_yri_to_eur373}

# -------------------- #
# EUR278 to YRI

source ${BASH_schedule_eur278_to_yri}

# -------------------- #
# EUR278 to FIN

source ${BASH_schedule_eur278_to_fin}


## -------------------- #
# now compile results
# will submit this as qsub job
# ${BASH_compile_results_continentalpop}

# binaries
Rscript=${Rscript}
qsub_variable_list="Rscript=${Rscript}"

# external scripts
R_compile_results_continentalpop=${R_compile_results_continentalpop}
qsub_variable_list="${qsub_variable_list},R_compile_results_continentalpop=${R_compile_results_continentalpop}"

# directories
logdir=${logdir}
datadir="${HOME}/gala_sage/rnaseq/geuvadis/rnaseq"
qsub_variable_list="${qsub_variable_list},logdir=${logdir}"

# file paths
eur373_to_eur373_path="${resultsdir}/eur373/geuvadis_elasticnet_eur373_predictions.txt"
eur373_to_afr_path="${resultsdir}/eur373/geuvadis_elasticnet_eur373_predictinto_yri89.txt"
eur278_to_eur278_path="${resultsdir}/eur278/geuvadis_elasticnet_eur278_predictions.txt"
eur278_to_fin_path="${resultsdir}/eur278/geuvadis_elasticnet_eur278_predictinto_fin95.txt"
eur278_to_afr_path="${resultsdir}/eur278/geuvadis_elasticnet_eur278_predictinto_yri89.txt"
afr_to_eur373_path="${resultsdir}/yri89/geuvadis_elasticnet_yri89_predictinto_eur373.txt"
afr_to_afr_path="${resultsdir}/yri89/geuvadis_elasticnet_yri89_predictions.txt"
fin_ids_path="${datadir}/geuvadis.fin.sampleids.txt"
eur278_ids_path="${datadir}/geuvadis.eur278.sampleids.txt"
output_predictions="${analysisdir}/geuvadis.continentalpop.predictions.txt"
output_r2="${analysisdir}/geuvadis.continentalpop.r2.txt"
poscorr="${analysisdir}/geuvadis.continentalpop.r2.poscorr.txt"
common_genes_poscorr="${analysisdir}/geuvadis.continentalpop.r2.poscorr.commongenes.txt"
qsub_variable_list="${qsub_variable_list},eur373_to_eur373_path=${eur373_to_eur373_path},eur373_to_afr_path=${eur373_to_afr_path},eur278_to_eur278_path=${eur278_to_eur278_path},eur278_to_fin_path=${eur278_to_fin_path},eur278_to_afr_path=${eur278_to_afr_path},afr_to_eur373_path=${afr_to_eur373_path},eur278_ids_path=${eur278_ids_path},fin_ids_path=${fin_ids_path},afr_to_afr_path=${afr_to_afr_path},output_results=${output_results},output_r2=${output_r2}"

# variables
h_rt="00:29:59"

# name of job?
compile_results="geuvadis.analyze.results"

# must hold this job against *all* other jobs previously scheduled
build_models="glmnet.${glmmethod}.eur373.to.yri89,glmnet.${glmmethod}.eur278.to.yri89,glmnet.${glmmethod}.eur278.to.fin95,glmnet.${glmmethod}.yri89.to.eur373"
collect_weights="glmnet.${glmmethod}.collect.weights.eur373.to.yri89,glmnet.${glmmethod}.collect.weights.eur278.to.yri89,glmnet.${glmmethod}.collect.weights.eur278.to.fin95,glmnet.${glmmethod}.collect.weights.yri89.to.eur373"
collect_weights="glmnet.${glmmethod}.collect.weights.eur373"
postprocess_results="glmnet.${glmmethod}.postprocess.eur373.to.yri89,glmnet.${glmmethod}.postprocess.eur278.to.yri89,glmnet.${glmmethod}.postprocess.eur278.to.fin95,glmnet.${glmmethod}.postprocess.yri89.to.eur373"
hold_jid="${build_models},${collect_weights},${postprocess_results}"

# execute
qsub -N ${compile_results} \
    -hold_jid ${hold_jid} \ 
    -v ${qsub_variable_list} \
    -o ${logdir} \
    -e ${logdir} \
    -l mem_free=${memory_limit} \
    -l h_rt=${h_rt} \
    ${BASH_compile_results_continentalpop}
