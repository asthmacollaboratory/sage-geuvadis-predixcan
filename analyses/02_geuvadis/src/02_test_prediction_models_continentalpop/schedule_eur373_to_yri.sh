#!/usr/bin/env bash

# set variables to population-specific parameters
# point prediction, weight, lambda, expression, subject files to EUR
# set alternative population parameters to YRI

h_rt="23:59:59"
pop="eur373"
subjectids=${subjectids_eur}
exprfile=${exprfile_eur}
phenofile=${phenofile_eur}
outdir=${outdir_eur}
resultssubdir=${resultssubdir_eur}
resultsdir=${resultsdir_eur}
predictionfile=${predictionfile_eur}
predictionfile_samepop=${predictionfile_eur2eur}
predictionfile_header=${predictionfile_header_eur373}
lambdafile=${lambdafile_eur}
weightsfile=${weightsfile_eur}
newweightsfile=${newweightsfile_eur}
out_lm_file=${out_lm_file_eur}
out_genelm_file=${out_genelm_file_eur}
num_pred_file=${num_pred_file_eur}
nsamples=${nsamples_eur}
h2file=${h2file_eur}
h2file_null=${h2file_null_eur}
nfolds=${nfolds_eur}
seed=${seed}


altpop="yri89"
altpop_exprfile=${exprfile_yri}
subjectids_altpop=${subjectids_yri}
predictionfile_altpop=${predictionfile_eur2yri}
altpop_out_lm_file=${out_lm_file_eur2yri}
altpop_out_genelm_file=${out_genelm_file_eur2yri}

predictionfile_samepop=${predictionfile_eur2eur}
samepop_out_lm_file=${out_lm_file_eur2eur}
samepop_out_genelm_file=${out_genelm_file_eur2eur}


# schedule jobs
source ${BASH_schedule_jobs}
