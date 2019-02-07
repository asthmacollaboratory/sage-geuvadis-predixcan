#!/usr/bin/env bash       # -- what is the language of this shell?
#                         # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash           # -- the shell for the job
#$ -r y                   # -- tell the system that if a job crashes, it should be restarted
#$ -j y                   # -- tell the system that the STDERR and STDOUT should be joined
#$ -l arch=linux-x64      # -- SGE resources (CPU type)
# ==========================================================================================
# copyright Asthma Collboratory (2018)
# coded by Kevin L. Keys
# ==========================================================================================


# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -o errexit # set -e, script will exit on error
set -o nounset # set -u, script will exit if it sees an uninitialized variable
set -o xtrace  # set -x, script will track which command is currently running
#ulimit -c 0 # user limits: -c covers the max size of core files created


# ==========================================================================================
# script variables (passed from QSUB command)
# ==========================================================================================

# binaries
Rscript=${Rscript}

# external scripts
R_postprocess_weights=${R_postprocess_weights}

# directories
logdir=${logdir}

# file paths
weightsfile=${weightsfile}
num_pred_file=${num_pred_file}
predictionfile_altpop=${predictionfile_altpop}
altpop_exprfile=${altpop_exprfile}
altpop_out_lm_file=${altpop_out_lm_file}
altpop_out_genelm_file=${altpop_out_genelm_file}
predictionfile=${predictionfile}
predictionfile_samepop=${predictionfile_samepop}
samepop_out_lm_file=${samepop_out_lm_file}
samepop_out_genelm_file=${samepop_out_genelm_file}
exprfile=${exprfile}
out_lm_file=${out_lm_file}
out_genelm_file=${out_genelm_file}

# other variables
glmmethod=${glmmethod}
pop=${pop}
nsamples=${nsamples}
discard_ratio=${discard_ratio}


# ==========================================================================================
# executable code
# ==========================================================================================

# start by noting current date, time, and process hostname
echo "Date: $(date)"
echo "Host name: $(hostname)"

# postprocess the weights file
$Rscript $R_postprocess_weights \
    --beta-file ${weightsfile} \
    --discard-ratio ${discard_ratio} \
    --num-predictions-file ${num_pred_file} \
    --num-samples ${nsamples} \
    --prediction-file ${predictionfile} \
    --expression-file ${exprfile} \
    --out-lm-file ${out_lm_file} \
    --out-genelm-file ${out_genelm_file} \
    --test-pop-prediction-file ${predictionfile_altpop} \
    --test-pop-expression-file ${altpop_exprfile} \
    --test-pop-out-lm-file ${altpop_out_lm_file} \
    --test-pop-out-genelm-file ${altpop_out_genelm_file} \
    --train-pop-prediction-file ${predictionfile_samepop} \
    --train-pop-out-lm-file ${samepop_out_lm_file} \
    --train-pop-out-genelm-file ${samepop_out_genelm_file}

# query return value of previous command
RETVAL=$?

# if return value is not 0, then previous command did not exit correctly
# create a status file notifying of error
# in contrary case, notify of success
if [ "${RETVAL}" -ne "0" ];
then
    echo "ERROR" > ${logdir}/status.${glmmethod}.postprocess.${pop}
else
    echo "SUCCESS" > ${logdir}/status.${glmmethod}.postprocess.${pop}
fi

# output end time
echo -e "End timestamp: $(date)"
