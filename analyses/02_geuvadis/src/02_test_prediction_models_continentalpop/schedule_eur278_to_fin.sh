#!/usr/bin/env bash

# set variables to population-specific parameters
# point prediction, weight, lambda, expression, subject files to EUR
# set alternative population parameters to FINN

h_rt="23:59:59"
pop="eur278"
subjectids=${subjectids_eur278}
exprfile=${exprfile_eur278}
phenofile=${phenofile_eur278}
outdir=${outdir_eur278}
resultssubdir=${resultssubdir_eur278}
resultsdir=${resultsdir_eur278}
predictionfile=${predictionfile_eur278}
predictionfile_samepop=${predictionfile_eur278toeur278}
predictionfile_header=${predictionfile_header_eur278}
lambdafile=${lambdafile_eur278}
weightsfile=${weightsfile_eur278}
newweightsfile=${newweightsfile_eur278}
out_lm_file=${out_lm_file_eur278}
out_genelm_file=${out_genelm_file_eur278}
num_pred_file=${num_pred_file_eur278}
nsamples=${nsamples_eur278}
h2file=${h2file_eur278}
h2file_null=${h2file_null_eur278}
nfolds=${nfolds_eur278}
seed=${seed}


altpop="fin95"
altpop_exprfile=${exprfile_fin}
subjectids_altpop=${subjectids_fin}
predictionfile_altpop=${predictionfile_eur278tofin}
altpop_out_lm_file=${out_lm_file_eur278tofin}
altpop_out_genelm_file=${out_genelm_file_eur278tofin}


# schedule jobs
source ${BASH_schedule_jobs}
