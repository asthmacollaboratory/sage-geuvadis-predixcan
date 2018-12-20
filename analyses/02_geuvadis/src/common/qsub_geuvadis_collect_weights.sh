#!/usr/bin/env bash       # -- what is the language of this shell
#$ -S /bin/bash           # -- the shell for the job
#$ -r y                   # -- tell the system that if a job crashes, it should be restarted
#$ -j y                   # -- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=1G         # -- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64      # -- SGE resources (CPU type)
#$ -l h_rt=0:10:00        # -- runtime limit in hours
# ==========================================================================================
# copyright Asthma Collboratory (2018)
# coded by Kevin L. Keys
#
# This script gathers results from a run of compute_new_predixcan_weights.sh.
# It produces two files:
# -- $resultsfile, which compiles prediction results for the current glmnet run;
# -- $weightsfile, which compiles nonzero prediction weights for the same glmnet run
# ==========================================================================================

# user limits: -c max size of core files created
echo "Start time: $(date)"
echo "Host name: $(hostname)"
ulimit -c 0

# script variables from qsub
weightsfile=${weightsfile}
glmnetdir=${glmnetdir}
glmmethod=${glmmethod}
logdir=${logdir}
lambdafile=${lambdafile}
predictionfile=${predictionfile}
outdir=${outdir}
resultsdir=${resultsdir}
resultssubdir=${resultssubdir}
tmpdir=${tmpdir}
pop=${pop}
altpop=${altpop}
predictionfile_altpop=${predictionfile_altpop}
predictionfile_samepop=${predictionfile_samepop}
predictionfile_header=${predictionfile_header}
h2file=${h2file}
h2file_null=${h2file_null}

# check variables for debugging
echo "Saving results to the following locations:"
echo -e "\tresultsfile = ${resultsfile}"
echo -e "\tweightsfile = ${weightsfile}"
echo -e "\tglmnetdir = ${glmnetdir}"
echo -e "\tglmmethod = ${glmmethod}"
echo -e "\tlogdir = ${logdir}"

# put header on prediction file
# then append results to file
rm -f ${predictionfile}
touch ${predictionfile}
echo -e "${predictionfile_header}" > ${predictionfile}
cat ${resultssubdir}/geuvadis_predictions_${glmmethod}_ENSG* | sort --temporary-directory ${tmpdir} >> ${predictionfile}

# check previous command
RETVAL=$?

# report on previous command
echo "exit status after making prediction file: ${RETVAL}"

# clean up temporary directory
rm -f ${tmpdir}/*

# put header on lambda file
# then computed weights to file
# ensure when concatenating results that we discard the individual file headers
rm -f ${lambdafile}
touch ${lambdafile}
echo -e "Gene\tHeld_out_Sample\tMean_MSE\tLambda\tPredicted_Expr" > ${lambdafile}
cat ${resultssubdir}/geuvadis_lambdas_${glmmethod}_ENSG* | grep -v "Gene" | sort --temporary-directory ${tmpdir} >> ${lambdafile}

# check previous command
# compound with penultimate one
let "RETVAL+=$?"
echo "exit status after making prediction,lambda files: ${RETVAL}"

# clean up temporary directory
rm -f ${tmpdir}/*

# now compile weights file
rm -f ${weightsfile}
touch ${weightsfile}
echo -e "Gene\tHeld_out_Sample\tSNP\tA1\tA2\tBeta" > ${weightsfile}
cat ${resultssubdir}/geuvadis_weights_${glmmethod}_ENSG* | grep -v "Gene" | fgrep -v "NA"$'\t'"NA" | sort --temporary-directory ${tmpdir} >> ${weightsfile}

# report on previous commands
let "RETVAL+=$?"
echo "exit status after making prediction,lambda,weights files: ${RETVAL}"

# compile external and internal prediction files
rm -f ${predictionfile_altpop}
touch ${predictionfile_altpop}
echo -e "Gene\tSubjectID\tPrediction_${pop}_into_${altpop}" > ${predictionfile_altpop}
cat ${resultssubdir}/geuvadis_predictinto_${altpop}_${glmmethod}_ENSG* | grep -F -v "SubjectID" | grep -F "ENS" | sort --temporary-directory ${tmpdir} >> ${predictionfile_altpop}

rm -f ${predictionfile_samepop}
touch ${predictionfile_samepop}
echo -e "Gene\tSubjectID\tPrediction_${pop}_into_${pop}" > ${predictionfile_samepop}
cat ${resultssubdir}/geuvadis_predictinto_${pop}_${glmmethod}_ENSG* | grep -F -v "SubjectID" | grep -F "ENS" | sort --temporary-directory ${tmpdir} >> ${predictionfile_samepop}

# finally, compile null and genotype heritabilities for all genes
rm -f ${h2file}
touch ${h2file}
echo -e "Gene\th2\th2_CI_lower\th2_CI_upper" > ${h2file}
cat ${resultssubdir}/geuvadis_h2_${pop}_ENS* | grep -F "ENS" | sort --temporary-directory ${tmpdir} | uniq >> ${h2file}

let "RETVAL+=$?"

rm -f ${h2file_null}
touch ${h2file_null}
echo -e "Gene\th2\tStdErr" > ${h2file_null}
cat ${resultssubdir}/geuvadis_h2_null_${pop}_ENS* | grep -F "ENS" | sort --temporary-directory ${tmpdir} | uniq >> ${h2file_null}


# report on previous commands
let "RETVAL+=$?"
echo "exit status after making prediction, lambda, weights, and external prediction files: ${RETVAL}"


# clean up temporary directory
rm -f ${tmpdir}/*

# now it is safe to clean up scratch if desired
#rm -rf $outdir

if [ "${RETVAL}" -ne "0" ];
then
    echo "ERROR" > ${logdir}/status.collectweights.${glmmethod}.${pop}
else
    echo "SUCCESS" > ${logdir}/status.collectweights.${glmmethod}.${pop}
fi

# announce stop time
echo "Stop time: $(date)"
