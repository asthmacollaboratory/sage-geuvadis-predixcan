#!/usr/bin/env bash
# ==============================================================================================================================
# copyright Asthma Collboratory (2018)
# coded by Kevin L. Keys
#
# This script contains all analysis-wide variables used in testing crosspopulation transcriptome imputation
# with GEUVADIS genotype-expression data.
#
# ==============================================================================================================================

# ==============================================================================================================================
# configure BASH
# ==============================================================================================================================
# these are used for debugging any undefined variables
#set -e
#set -u


# ==============================================================================================================================
# directories
# ==============================================================================================================================

# script static directories
thisdir="$(dirname $(readlink -f $0))"

analysisdir="${thisdir}/../../analysis"
datadir="${analysisdir}/data"
datafiles_dir="${thisdir}/../../datafiles"
rnaseqdir="${datadir}/rnaseq"
genodir="${datadir}/genotypes/plink"
gctadir="${datadir}/gcta"
resultsdir="${analysisdir}/results"
logdir="${resultsdir}/log"
commondir="${thisdir}/../common"

# directories that follow rely on glmnet settings
glmmethod="elasticnet"
alpha="0.5"

# on UCSF QB3, the scratch directory is a fixed path
# change this to wherever SGE results can be stored temporarily
# it should have plenty of space, say >100Gb
scratchdir="/scratch/kkeys"

outdir_eur="${scratchdir}/${glmmethod}/genes/eur373"
outdir_eur278="${scratchdir}/${glmmethod}/genes/eur278"
outdir_yri="${scratchdir}/${glmmethod}/genes/yri89"
outdir_tsi="${scratchdir}/${glmmethod}/genes/tsi"
outdir_gbr="${scratchdir}/${glmmethod}/genes/gbr"
outdir_fin="${scratchdir}/${glmmethod}/genes/fin"
outdir_ceu="${scratchdir}/${glmmethod}/genes/ceu"

resultsdir_eur="${resultsdir}/${glmmethod}/eur373"
resultsdir_eur278="${resultsdir}/${glmmethod}/eur278"
resultsdir_yri="${resultsdir}/${glmmethod}/yri89"
resultsdir_ceu="${resultsdir}/${glmmethod}/ceu92"
resultsdir_gbr="${resultsdir}/${glmmethod}/gbr96"
resultsdir_fin="${resultsdir}/${glmmethod}/fin95"
resultsdir_tsi="${resultsdir}/${glmmethod}/tsi93"

resultsdir_crosspop="${resultsdir}/${glmmethod}/crosspop"

resultsdir_ceu89="${resultsdir_crosspop}/ceu89"
resultsdir_tsi89="${resultsdir_crosspop}/tsi89"
resultsdir_gbr89="${resultsdir_crosspop}/gbr89"
resultsdir_fin89="${resultsdir_crosspop}/fin89"

resultssubdir_eur="${resultsdir_eur}/results"
resultssubdir_eur278="${resultsdir_eur278}/results"
resultssubdir_yri="${resultsdir_yri}/results"
resultssubdir_ceu="${resultsdir_ceu}/results"
resultssubdir_gbr="${resultsdir_gbr}/results"
resultssubdir_fin="${resultsdir_fin}/results"
resultssubdir_tsi="${resultsdir_tsi}/results"
resultssubdir_ceu89="${resultsdir_ceu89}/results"
resultssubdir_gbr89="${resultsdir_gbr89}/results"
resultssubdir_fin89="${resultsdir_fin89}/results"
resultssubdir_tsii89="${resultsdir_tsi89}/results"

tmpdir="${analysisdir}/tmp"

crosspop_dir="${rnaseqdir}/crosspop"

# make output and results directories, in case they doesn't exist
# do this before defining filepaths and paths to binaries
mkdir -p ${outdir_eur}
mkdir -p ${outdir_eur278}
mkdir -p ${outdir_yri}
mkdir -p ${resultsdir_eur}
mkdir -p ${resultsdir_eur278}
mkdir -p ${resultsdir_yri}
mkdir -p ${resultsdir_ceu}
mkdir -p ${resultsdir_gbr}
mkdir -p ${resultsdir_fin}
mkdir -p ${resultsdir_tsi}
mkdir -p ${resultssubdir_eur}
mkdir -p ${resultssubdir_eur278}
mkdir -p ${resultssubdir_yri}
mkdir -p ${resultssubdir_ceu}
mkdir -p ${resultssubdir_gbr}
mkdir -p ${resultssubdir_fin}
mkdir -p ${resultssubdir_tsi}
mkdir -p ${tmpdir}
mkdir -p ${crosspop_dir}


# ==============================================================================================================================
# filepaths
# ==============================================================================================================================

bedfile_pfx="${genodir}/GEUVADIS.ALLCHR.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.rsq_0.8_maf_0.01_hwe_0.00001_geno_0.05"

exprfile_eur="${rnaseqdir}/geuvadis.eur373.RPKM.invnorm.txt"
exprfile_eur278="${rnaseqdir}/geuvadis.eur278.RPKM.invnorm.txt"
exprfile_yri="${rnaseqdir}/geuvadis.yri89.RPKM.invnorm.txt"
exprfile_fin="${rnaseqdir}/geuvadis.fin95.RPKM.invnorm.txt"

phenofile_eur="${rnaseqdir}/geuvadis.eur373.RPKM.invnorm.pheno"
phenofile_eur278="${rnaseqdir}/geuvadis.eur278.RPKM.invnorm.pheno"
phenofile_yri="${rnaseqdir}/geuvadis.yri89.RPKM.invnorm.pheno"
phenofile_fin="${rnaseqdir}/geuvadis.fin95.RPKM.invnorm.pheno"

genelist="${datafiles_dir}/human_ens_GRCh37_genechrpos_plusmin500kb.txt"

subjectids_eur="${rnaseqdir}/geuvadis.eur373.sampleids.txt"
subjectids_yri="${rnaseqdir}/geuvadis.yri89.sampleids.txt"
subjectids_fin="${rnaseqdir}/geuvadis.fin95.sampleids.txt"
subjectids_eur278="${rnaseqdir}/geuvadis.eur278.sampleids.txt"

predictionfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictions.txt"
predictionfile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictions.txt"
predictionfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictions.txt"

predictionfile_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89.txt"
predictionfile_eur2eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_eur373.txt"
predictionfile_eur278toyri="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_yri89.txt"
predictionfile_eur278tofin="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_fin95.txt"
predictionfile_eur278toeur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_eur278.txt"
predictionfile_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373.txt"
predictionfile_yri2yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_yri89.txt"

lambdafile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_lambdas.txt"
lambdafile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_lambdas.txt"
lambdafile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_lambdas.txt"

weightsfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_weights.txt"
weightsfile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_weights.txt"
weightsfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_weights.txt"

num_pred_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_numpred.txt"
num_pred_file_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_numpred.txt"
num_pred_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_numpred.txt"

newweightsfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_weights_noNA.txt"
newweightsfile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_weights_noNA.txt"
newweightsfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_weights_noNA.txt"

out_lm_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_lm_predvmeas_results.txt"
out_lm_file_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_lm_predvmeas_results.txt"
out_lm_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_lm_predvmeas_results.txt"

out_genelm_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_genelm_predvmeas_results.txt"
out_genelm_file_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_genelm_predvmeas_results.txt"
out_genelm_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_genelm_predvmeas_results.txt"

out_lm_file_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89_lm_predvmeas_results.txt"
out_lm_file_eur2eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_eur373_lm_predvmeas_results.txt"
out_lm_file_eur278toeur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_eur278_lm_predvmeas_results.txt"
out_lm_file_eur278tofin="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_fin95_lm_predvmeas_results.txt"
out_lm_file_eur278toyri="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_yri89_lm_predvmeas_results.txt"
out_lm_file_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373_lm_predvmeas_results.txt"
out_lm_file_yri2yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_yri89_lm_predvmeas_results.txt"

out_genelm_file_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89_genelm_predvmeas_results.txt"
out_genelm_file_eur2eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_eur373_genelm_predvmeas_results.txt"
out_genelm_file_eur278toeur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_eur278_genelm_predvmeas_results.txt"
out_genelm_file_eur278toyri="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_yri89_genelm_predvmeas_results.txt"
out_genelm_file_eur278tofin="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_fin95_genelm_predvmeas_results.txt"
out_genelm_file_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373_genelm_predvmeas_results.txt"
out_genelm_file_yri2yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_yri89_genelm_predvmeas_results.txt"

h2file_eur="${resultsdir_eur}/geuvadis_h2_eur373.txt"
h2file_null_eur="${resultsdir_eur}/geuvadis_h2_null_eur373.txt"
h2file_eur278="${resultsdir_eur278}/geuvadis_h2_eur278.txt"
h2file_null_eur278="${resultsdir_eur278}/geuvadis_h2_null_eur278.txt"
h2file_yri="${resultsdir_yri}/geuvadis_h2_yri89.txt"
h2file_null_yri="${resultsdir_yri}/geuvadis_h2_null_yri89.txt"

# external script locations
R_compute_new_weights="${commondir}/geuvadis_compute_new_weights.R"
R_postprocess_weights="${commondir}/geuvadis_postprocess_weights.R"
R_predict_new_pop="${commondir}/geuvadis_predict_in_altpop.R"
R_subsample_pop="${commondir}/geuvadis_subsample_onepop.R"

BASH_schedule_jobs="${commondir}/geuvadis_qsub_jobs.sh"
BASH_compute_weights="${commondir}/qsub_geuvadis_compute_weights.sh"
BASH_collect_weights="${commondir}/qsub_geuvadis_collect_weights.sh"
BASH_postprocess_weights="${commondir}/qsub_geuvadis_postprocess_weights.sh"

# assumes that FIESTA is cloned into $HOME/git
# edit this as needed
PYTHON_fiesta="${HOME}/Git/albi/fiesta.py"

# binaries
#PLINK=$(whereis plink | awk '{print $2}')
PLINK="${HOME}/bin/plink" ## QB3 doesn't have system-wide PLINK
Rscript=$(whereis Rscript | awk '{print $2}')
GCTA=$(whereis gcta64| awk '{print $2}')
PYTHON=$(whereis python2.7| awk '{print $2}')

# script variables
maf="0.01"
hwe="0.0001"
geno="0.03"
nthreads=1
memory_limit="5G"
memory_limit_mb="5000" # manually coordinate this with $memory_limit!!!
scratch_memory="2G"
h_rt="23:59:59"
discard_ratio="0.5" # desired min ratio of samples with nonmissing LOOCV predictions, used in postprocessing, e.g. 0.5 = 50%
nsamples_eur="373"
nsamples_eur278="278"
nsamples_yri="89"
nfolds_eur="10"
nfolds_eur278="10"
nfolds_yri="89"
seed=2018

# need headers for prediction files
predictionfile_header_eur373="$(head -n 1 ${exprfile_eur})"
predictionfile_header_yri="$(head -n 1 ${exprfile_yri})"
predictionfile_header_eur278="$(head -n 1 ${exprfile_eur278})"

# ==============================================================================================================================
# start script
# ==============================================================================================================================

# -------------------- #
echo "Start time: $(date)"

# -------------------- #
# make necessary output directories
if [[ ! -d "${logdir}" ]]; then mkdir -p ${logdir}; fi


# how many genes are in this list?
# remember that this file has a header so line count is off by 1
# another way to get correct $nGenes is to follow wc with
# > let "nGenes=nGenes-1"
nGenes=$(wc -l ${genelist} | cut -f 1 -d " ")

# uncomment next line for testing and debugging the pipeline
#nGenes=10
