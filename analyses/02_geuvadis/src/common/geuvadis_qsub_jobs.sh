#!/usr/bin/env bash
# ==============================================================================================================================
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys
#
# This script is called from the main GEUVADIS prediction model training pipeline. IT DOES NOT WORK OUTSIDE OF THE SCRIPT.
# The requisite variables in the following qsub arguments are set in the pipeline. DO NOT CHANGE THEM HERE.
#
# This script schedules three jobs:
#     1. compute new weights from a training set of genotypes and gene expression measures
#     2. collect the individual gene files and compile them into a handful of files
#     3. postprocess the results, generate imputation R2 and correlations
# ==============================================================================================================================

# ==============================================================================================================================
# 0. variables for this script only 
# ==============================================================================================================================

# job names
build_models="glmnet.${glmmethod}.${pop}"
collect_weights="glmnet.${glmmethod}.collect.weights.${pop}"
postprocess_results="glmnet.${glmmethod}.postprocess.${pop}"

# ==============================================================================================================================
# 1. compute new weights from a training set of genotypes and gene expression measures
# ==============================================================================================================================

# build list of variables to pass to QSUB option -v
# will do this once for each of the three scripts below

# start with binaries
Rscript=${Rscript}
FIESTA=${PYTHON_fiesta}
PYTHON=${PYTHON}
GCTA=${GCTA}
PLINK=${PLINK}
qsub_variable_list="Rscript=${Rscript},FIESTA=${PYTHON_fiesta},PYTHON=${PYTHON},GCTA=${GCTA},PLINK=${PLINK}"

# add external scripts
R_predict_new_pop=$R_predict_new_pop
R_compute_new_weights=$R_compute_new_weights
qsub_variable_list="${qsub_variable_list},R_predict_new_pop=$R_predict_new_pop,R_compute_new_weights=$R_compute_new_weights"

# add directories
logdir=${logdir}
outdir=${outdir}
resultssubdir=${resultssubdir}
resultsdir=${resultsdir}
gctadir=${gctadir}
imputegenodir=${imputegenodir}
tmpdir=${tmpdir}
qsub_variable_list="${qsub_variable_list},logdir=${logdir},outdir=${outdir},resultssubdir=${resultssubdir},resultsdir=${resultsdir},gctadir=${gctadir},imputegenodir=${imputegenodir},tmpdir=${tmpdir}"

# add filepaths
genelist=${genelist}
subjectids=${subjectids}
exprfile=${exprfile}
subjectids_altpop=${subjectids_altpop}
phenofile=${phenofile}
qsub_variable_list="${qsub_variable_list},genelist=${genelist},subjectids=${subjectids},exprfile=${exprfile},subjectids_altpop=${subjectids_altpop},phenofile=${phenofile}"

# add numeric or string variables
alpha=${alpha}
glmmethod=${glmmethod}
maf=${maf}
hwe=${hwe}
nthreads=${nthreads}
memory_limit_mb=${memory_limit_mb}
pop=${pop}
altpop=${altpop}
nfolds=${nfolds}
qsub_variable_list="${qsub_variable_list},alpha=${alpha},glmmethod=${glmmethod},maf=${maf},hwe=${hwe},nthreads=${nthreads},memory_limit_mb=${memory_limit_mb},pop=${pop},altpop=${altpop},nfolds=${nfolds}"

# execute
qsub -N ${build_models} \
     -v ${qsub_variable_list} \
     -t 1-${nGenes} \
     -e ${logdir} \
     -o ${logdir} \
     -l mem_free=${memory_limit} \
     -l scratch=${scratch_memory} \
     -l h_rt=${h_rt} \
     ${BASH_compute_weights}


# ==============================================================================================================================
# 2. collect the individual gene files and compile them into a handful of files
# ==============================================================================================================================

# beefy calculations done, may not need as much compute time
# can set $h_rt to < 30 min to schedule in QB3 speedy queue
# but nota bene: jobs in short.q killed at 30 min regardless of completion status!!!
#h_rt="00:29:59"
h_rt="06:00:00"

# build variable list
# start with directories
glmnetdir=${glmnetdir}
logdir=${logdir}
outdir=${outdir}
resultsdir=${resultsdir}
resultssubdir=${resultssubdir}
tmpdir=${tmpdir}
qsub_variable_list="glmnetdir=${glmnetdir},logdir=${logdir},outdir=${outdir},resultsdir=${resultsdir},resultssubdir=${resultssubdir},tmpdir=${tmpdir}"

# add file paths}
weightsfile=${weightsfile}
predictionfile=${predictionfile}
lambdafile=${lambdafile}
predictionfile_altpop=${predictionfile_altpop}
predictionfile_samepop=${predictionfile_samepop}
h2file=${h2file}
h2file_null=${h2file_null}
qsub_variable_list="${qsub_variable_list},weightsfile=${weightsfile},predictionfile=${predictionfile},lambdafile=${lambdafile},predictionfile_altpop=${predictionfile_altpop},predictionfile_samepop=${predictionfile_samepop},h2file=${h2file},h2file_null=${h2file_null}"

# add remaining variables
glmmethod=${glmmethod}
pop=${pop}
altpop=${altpop}
predictionfile_header=${predictionfile_header}
qsub_variable_list="${qsub_variable_list},glmmethod=${glmmethod},pop=${pop},altpop=${altpop},predictionfile_header=${predictionfile_header}"

# execute
qsub -N ${collect_weights} \
    -hold_jid ${build_models} \ 
    -v ${qsub_variable_list} \
    -o ${logdir} \
    -e ${logdir} \
    -l mem_free=${memory_limit} \
    -l h_rt=${h_rt} \
    ${BASH_collect_weights}

# ==============================================================================================================================
# 3. postprocess the results, generate imputation R2 and correlations
# ==============================================================================================================================

# same deal as before

# binaries
Rscript=${Rscript}
qsub_variable_list="Rscript=${Rscript}"

# external scripts
R_glmnet_postprocess=${R_glmnet_postprocess}
qsub_variable_list="${qsub_variable_list},R_glmnet_postprocess=${R_glmnet_postprocess}"

# directories
glmnetdir=${glmnetdir}
logdir=${logdir}
qsub_variable_list="${qsub_variable_list},glmnetdir=${glmnetdir},logdir=${logdir}"

# file paths
weightsfile=${weightsfile}
newweightsfile=${newweightsfile}
num_pred_file=${num_pred_file}
predictionfile_altpop=${predictionfile_altpop}
altpop_exprfile=${altpop_exprfile}
altpop_out_lm_file=${altpop_out_lm_file}
altpop_out_genelm_file=${altpop_out_genelm_file}
predictionfile=${predictionfile}
exprfile=${exprfile}
out_lm_file=${out_lm_file}
out_genelm_file=${out_genelm_file}
qsub_variable_list="${qsub_variable_list},weightsfile=${weightsfile},newweightsfile=${newweightsfile},num_pred_file=${num_pred_file},predictionfile_altpop=${predictionfile_altpop},altpop_exprfile=${altpop_exprfile},altpop_out_lm_file=${altpop_out_lm_file},altpop_out_genelm_file=${altpop_out_genelm_file},predictionfile=${predictionfile},exprfile=${exprfile},out_lm_file=${out_lm_file},out_genelm_file=${out_genelm_file}"

# variables
glmmethod=${glmmethod}
pop=${pop}
nsamples=${nsamples}
discard_ratio=${discard_ratio}
qsub_variable_list="${qsub_variable_list},glmmethod=${glmmethod},pop=${pop},nsamples=${nsamples},discard_ratio=${discard_ratio}"

# execute
qsub -N ${postprocess_results} \
    -hold_jid ${build_models},${collect_weights} \ 
    -v ${qsub_variable_list} \ 
    -o ${logdir} \
    -e ${logdir} \
    -l mem_free=${memory_limit} \
    -l h_rt=${h_rt} \
    ${BASH_postprocess_weights}
