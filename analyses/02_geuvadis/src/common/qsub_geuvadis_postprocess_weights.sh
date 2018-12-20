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

# script variables
glmmethod=${glmmethod}
weightsfile=${weightsfile}
glmnetdir=${glmnetdir}
#newweightsfile=${newweightsfile}
Rscript=${Rscript}
R_glmnet_postprocess=${R_glmnet_postprocess}
discard_ratio=${discard_ratio}
num_pred_file=${num_pred_file}
nsamples=${nsamples}
predictionfile=${predictionfile}
exprfile=${exprfile}
out_lm_file=${out_lm_file}
out_genelm_file=${out_genelm_file}
logdir=${logdir}
pop=${pop}
predictionfile_altpop=${predictionfile_altpop}
altpop_exprfile=${altpop_exprfile}
altpop_out_lm_file=${altpop_out_lm_file}
altpop_out_genelm_file=${altpop_out_genelm_file}

# output start time
echo -e "Start timestamp: $(date)"

# postprocess the weights file
#$Rscript $R_glmnet_postprocess $weightsfile $newweightsfile $discard_ratio $num_pred_file $nsamples $predictionfile $exprfile $out_lm_file $out_genelm_file $predictionfile_altpop $altpop_exprfile $altpop_out_lm_file $altpop_out_genelm_file
$Rscript $R_glmnet_postprocess \
    --beta-file ${weightsfile} \
    --discard-ratio ${$discard_ratio} \
    --num-predictions-file ${num_pred_file} \
    --num_samples ${nsamples} \
    --prediction-file ${predictionfile} \
    --expression-file ${exprfile} \
    --out-lm-file ${out_lm_file} \
    --out-genelm-file ${out_genelm_file} \
    --test-pop-prediction-file ${predictionfile_altpop} \
    --test-pop-expression-file ${altpop_exprfile} \
    --test-pop-out-lm-file ${altpop_out_lm_file} \
    --test-pop-out-genelm-file ${altpop_out_genelm_file}

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
