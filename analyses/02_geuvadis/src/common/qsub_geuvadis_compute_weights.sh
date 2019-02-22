#!/usr/bin/env bash                # -- what is the language of this shell?
#                                  # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    # -- the shell for the job
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined
#$ -l arch=linux-x64               # -- SGE resources (CPU type)
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
# This script prepares and executes training and testing of prediction models using
# GEUVADIS gene expression data.
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
PLINK=${PLINK}
Rscript=${Rscript}

# external scripts
R_compute_new_weights=${R_compute_new_weights}
R_predict_new_pop=${R_predict_new_pop}

# directories
logdir=${logdir}
outdir=${outdir}
gctadir=${gctadir}
resultsdir=${resultsdir}
resultssubdir=${resultssubdir}
tmpdir=${tmpdir}

# file paths
genelist=${genelist}
subjectids=${subjectids}
exprfile=${exprfile}
subjectids_altpop=${subjectids_altpop}
bedfile_pfx=${bedfile_pfx}

# other variables
alpha=${alpha}
glmmethod=${glmmethod}
maf=${maf}
hwe=${hwe}
nthreads=${nthreads}
memory_limit_mb=${memory_limit_mb}
pop=${pop}
altpop=${altpop}
nfolds=${nfolds}
seed=${seed}

# ==========================================================================================
# executable code
# ==========================================================================================

# start by noting current date, time, and process hostname
echo "Date: $(date)"
echo "Host name: $(hostname)"

# parse current gene
# NOTA BENE: in general BASH arrays are 0-indexed while SGE tasks are 1-indexed
# since $genelist lacks a header then $genelist is essentially 0-indexed
# must subtract 1 from $SGE_TASK_ID to match correct gene
# in general, be mindful when indexing BASH arrays with SGE task IDs
i=$(expr ${SGE_TASK_ID} - 1) || true  ## guard against exit status 0
read -a genes <<< $(cat ${genelist} | cut -d " " -f 1)
gene=${genes[${i}]}

echo "Preparing analysis of ${gene} in population ${pop}..."

# create directory path to gene folder
# make gene directory in case it doesn't exist
genepath="${outdir}/${gene}"
genopfx="${genepath}/${gene}"
mkdir -p ${genepath}

# create file paths to PLINK output
rawpath="${genopfx}.raw"
bimfile="${genopfx}.bim"

# also create paths for
genopfx_altpop="${genepath}/${gene}_${altpop}"
rawpath_altpop="${genopfx_altpop}.raw"

# create paths to output files for results
predictionfile="${resultssubdir}/geuvadis_predictions_${glmmethod}_${gene}.txt"
lambdafile="${resultssubdir}/geuvadis_lambdas_${glmmethod}_${gene}.txt"
weightsfile="${resultssubdir}/geuvadis_weights_${glmmethod}_${gene}.txt"
predictionfile_altpop="${resultssubdir}/geuvadis_predictinto_${altpop}_${glmmethod}_${gene}.txt"
predictionfile_samepop="${resultssubdir}/geuvadis_predictinto_${pop}_${glmmethod}_${gene}.txt"

# make note of output file paths, useful for debugging
echo "Will save output to the following files:"
echo -e "\tweightsfile = ${weightsfile}"
echo -e "\tpredictionfile = ${predictionfile}"
echo -e "\tlambdafile = ${lambdafile}"
echo -e "\tpredictionfile_altpop = ${predictionfile_altpop}"

# parse info for current gene
mygeneinfo=$(grep ${gene} ${genelist})
chr=$(echo ${mygeneinfo} | cut -f 2 -d " ")
startpos=$(echo ${mygeneinfo} | cut -f 3 -d " ")
endpos=$(echo ${mygeneinfo} | cut -f 4 -d " ")

echo "Information for current gene:"
echo -e "\tname: ${gene}"
echo -e "\tchr: ${chr}"
echo -e "\tstart: ${startpos}"
echo -e "\tend: ${endpos}"

# create a PLINK RAW file
# this codes the dosage format required for glmnet
# here we use the genome-wide GEUVADIS genotype data with rsIDs
echo "Subsetting genotypes in ${gene} for training pop ${pop}..."
$PLINK \
    --bfile ${bedfile_pfx} \
    --chr ${chr} \
    --from-bp ${startpos} \
    --to-bp ${endpos} \
    --maf ${maf} \
    --hwe ${hwe} \
    --recode A \
    --make-bed \
    --out ${genopfx} \
    --threads ${nthreads} \
    --memory ${memory_limit_mb} \
    --keep ${subjectids} #\
    #--silent ## turn this off first when debugging

# call glmnet script
# the method used depends on the alpha value:
# alpha = "0.5" --> elastic net regression
# alpha = "1.0" --> LASSO regression
# alpha = "0.0" --> ridge regression
echo "starting R script to compute new prediction weights..."
$Rscript $R_compute_new_weights \
    --genotype-dosage-file ${rawpath} \
    --expression-file ${exprfile} \
    --gene-name ${gene} \
    --prediction-output ${predictionfile} \
    --lambda-output ${lambdafile} \
    --alpha ${alpha} \
    --beta-output ${weightsfile} \
    --BIM-file ${bimfile} \
    --num-folds ${nfolds} \
    --random-seed ${seed}

# query return value of previous command
# will use this later to determine successful execution
RETVAL=$?

# get list of SNPs to subset in alternate population
echo "Constructing list of SNPs to extract in testing population ${altpop}..."
snps_to_extract="${tmpdir}/snps_to_extract_${altpop}_${gene}.txt"
cat ${weightsfile} | cut -f 3 | grep "rs" | sort | uniq > $snps_to_extract


# create a PLINK RAW file, but this time for the testing population
# this codes the dosage format required for glmnet
# here we use the genome-wide GEUVADIS genotype data with rsIDs
echo "Subsetting genotypes in ${gene} for the testing population ${altpop}..."
$PLINK \
    --bfile ${bedfile_pfx} \
    --recode A \
    --make-bed \
    --out ${genopfx_altpop} \
    --threads ${nthreads} \
    --memory ${memory_limit_mb} \
    --keep ${subjectids_altpop} \
    --extract ${snps_to_extract} #\
    #--silent ## turn this off first when debugging

# note the following commented PLINK options: why are they not used?
# we want to use all possible SNPs from training pop
# but genetic variation may not match that of testing pop
# tricky to filter on same MAF/HWE thresholds as a result since we may lose SNPs
# to be safe, just include all possible SNPs
#    --maf ${maf} \
#    --hwe ${hwe} \
#    --chr ${chr} \
#    --from-bp ${startpos} \
#    --to-bp ${endpos} \

# predict from training pop to testing pop
echo "Predicting into the testing population ${altpop}..."
$Rscript $R_predict_new_pop \
    --beta-file ${weightsfile} \
    --genotype-dosage-file ${rawpath_altpop} \
    --prediction-output ${predictionfile_altpop} \
    --gene-name ${gene}

# query return value of previous command
let "RETVAL+=$?" || true  ## need " || true" to satisfy "set -e"

# also perform prediction from training pop into itself
# different from true out-of-sample populations but systematically same as before
# want this to compare against quality of out-of-sample pops
echo "Predicting into the training population..."
$Rscript $R_predict_new_pop \
    --beta-file ${weightsfile} \
    --genotype-dosage-file ${rawpath} \
    --prediction-output ${predictionfile_samepop} \
    --gene-name ${gene}

# query return value of previous command
let "RETVAL+=$?" || true

# if return value is not 0, then previous command did not exit correctly
# create a status file notifying of error
# in contrary case, notify of success
if [ "$RETVAL" -ne "0" ];
then
  echo "ERROR" > ${logdir}/status.${glmmethod}.${gene}.${pop}
else
  echo "SUCCESS" > ${logdir}/status.${glmmethod}.${gene}.${pop}
fi

# tell job to report on itself
# also makes a useful handle for inspecting completed jobs
echo "job ${JOB_ID} report:"
qstat -j ${JOB_ID}
