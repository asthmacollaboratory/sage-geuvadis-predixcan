#!/usr/bin/env bash

# set variables to population-specific parameters
# point prediction, weight, lambda, expression, subject files to YRI
# set alternative population parameters to EUR

h_rt="12:00:00"  # yri has smaller sample size than eur373 --> less compute time necessary
pop="yri89"
subjectids=${subjectids_yri}
exprfile=${exprfile_yri}
phenofile=${phenofile_yri}
outdir=${outdir_yri}
resultssubdir=${resultssubdir_yri}
resultsdir=${resultsdir_yri}
predictionfile=${predictionfile_yri}
predictionfile_samepop=${predictionfile_yri2yri}
predictionfile_header=${predictionfile_header_yri}
lambdafile=${lambdafile_yri}
weightsfile=${weightsfile_yri}
newweightsfile=${newweightsfile_yri}
out_lm_file=${out_lm_file_yri}
out_genelm_file=${out_genelm_file_yri}
num_pred_file=${num_pred_file_yri}
nsamples=${nsamples_yri}
h2file=${h2file_yri}
h2file_null=${h2file_null_yri}
nfolds=${nfolds_yri}
seed=${seed}


altpop="eur373"
altpop_exprfile=${exprfile_eur}
subjectids_altpop=${subjectids_eur}
predictionfile_altpop=${predictionfile_yri2eur}
altpop_out_lm_file=${out_lm_file_yri2eur}
altpop_out_genelm_file=${out_genelm_file_yri2eur}


# schedule jobs
source ${BASH_schedule_jobs}
