#!/usr/bin/env bash
# ==========================================================================================
# copyright Asthma Collboratory (2018)
# coded by Kevin L. Keys
#
# This script computes PrediXcan weights from gEUVADIS transcriptome data.
#
# Call:
#
# ./compute_new_predixcan_weights_geuvadis.sh
# ==========================================================================================


# ==========================================================================================
# BASH script settings
# ==========================================================================================

set -o errexit # set -e, script will exit on error
set -o nounset # set -u, script will exit if it sees an uninitialized variable
#set -o xtrace  # set -x, script will track which command is currently running

# ==========================================================================================
# source variables
# ==========================================================================================

thisdir="$(dirname $(readlink -f $0))"
BASH_define_variables="${thisdir}/../common/geuvadis_variables.sh"

source $BASH_define_variables

# ==========================================================================================
# script variables
# ==========================================================================================

subsample_size="${nsamples_yri}"


# ==========================================================================================
# cross-population training
# ==========================================================================================

# in this section we will train from CEU, GBR, TSI, FIN, YRI into all other populations
# for example, we train in CEU and predict into other four populations
#
# key is that we downsample to smallest population. sample sizes:
# -- YRI = 89
# -- CEU = 92
# -- FIN = 95
# -- GBR = 96
# -- TSI = 93
#
# thus, we must downsample to 89 samples in each training population

# arrays of pops + their respective sample sizes
pops=("ceu" "tsi" "gbr" "fin" "yri")
popsizes=("92" "93" "96" "95" "89")

# set runtime to relatively high value
h_rt="23:59:59"

# loop over pops
for i in $(seq 0 4); do
#for i in $(seq 0 0); do

    # set variables for current pop
    # train in "pop", test in "notpop"
    pop="${pops[$i]}"  ## note absence of sample size; will add subsample size 89 to this later
    notpop="not${pop}"
    popsize="${popsizes[$i]}"
    altpop="${notpop}"
    subjectids="${datafiles_dir}/geuvadis.${pop}${popsize}.sampleids.txt"

    # first subsample the population in question
    $Rscript $R_subsample_pop \
        --population-name ${pop} \
        --EUR-RNA-file ${exprfile_eur} \
        --AFR-RNA-file ${exprfile_yri} \
        --sample-file ${subjectids} \
        --output-directory ${crosspop_dir} \
        --seed ${seed} \
        --subsample-size ${subsample_size}

    # rename population with correct subsampled population size
    pop="${pop}${subsample_size}"

    outdir="${scratchdir}/${glmmethod}/genes/${pop}"
    resultsdir="${resultsdir_crosspop}/${pop}"
    resultssubdir="${resultsdir}/results"

    mkdir -p ${outdir}
    mkdir -p ${resultsdir}
    mkdir -p ${resultssubdir}

    subjectids="${crosspop_dir}/geuvadis.${pop}.sampleids.txt"
    exprfile="${crosspop_dir}/geuvadis.${pop}.RPKM.invnorm.txt"
    predictionfile="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictions.txt"
    lambdafile="${resultsdir}/geuvadis_${glmmethod}_${pop}_lambdas.txt"
    weightsfile="${resultsdir}/geuvadis_${glmmethod}_${pop}_weights.txt"
    predictionfile_samepop="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${pop}.txt"
    samepop_out_lm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${pop}_lm_predvmeas_results.txt"
    samepop_out_genelm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${pop}_genelm_predvmeas_results.txt"
    num_pred_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_numpred.txt"
    nsamples=${subsample_size}
    out_lm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_lm_predvmeas_results.txt"
    newweightsfile="${resultsdir}/geuvadis_${glmmethod}_${pop}_weights_noNA.txt"
    out_genelm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_genelm_predvmeas_results.txt"
    h2file="${resultsdir}/geuvadis_h2_${pop}.txt"
    h2file_null="${resultsdir}/geuvadis_h2_null_${pop}.txt"
    nfolds=${nfolds_yri}
    phenofile="${crosspop_dir}/geuvadis.${pop}.RPKM.invnorm.pheno"

    altpop_exprfile="${crosspop_dir}/geuvadis.${notpop}.RPKM.invnorm.txt"
    subjectids_altpop="${crosspop_dir}/geuvadis.${notpop}.sampleids.txt"
    predictionfile_altpop="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${notpop}.txt"
    altpop_out_lm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${notpop}_lm_predvmeas_results.txt"
    altpop_out_genelm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${notpop}_genelm_predvmeas_results.txt"

    predictionfile_header="$(head -n 1 ${exprfile} | sed -e 's/,/\t/g')"


    source $BASH_schedule_jobs
done

# end script
echo "End time: $(date)"
