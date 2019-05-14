#!/usr/bin/env bash
# ==========================================================================================
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys and Donglei Hu
# ==========================================================================================

# ==========================================================================================
# script variables 
# ==========================================================================================

# directories
thisdir="$(dirname $(readlink -f $0))"

predixcan_sfx="predixcan.genotypefile.SAGE_hg19-1kgP3v5_panel.mixedpop_phased.eagle_39.pass.x-lowComplex.dp10-gq20.SNP.id.txt" # coordinate with Cesar output

analysisdir="${thisdir}/../../analysis"
datafiles_dir="${thisdir}/../../datafiles"

datadir="${analysisdir}/data"
dbdir="${analysisdir}/databases"

rnaseqdir="${datadir}/rnaseq"
genodir="${datadir}/genotypes"

resultsdir="${analysisdir}/results"
outdir="${resultsdir}"
predixcanoutdir="${resultsdir}/predixcan"

dgndir="${predixcanoutdir}/DGN"
gtexdir="${predixcanoutdir}/GTEx"
mesadir="${predixcanoutdir}/MESA"
logdir="${predixcanoutdir}/log"
predixcan_genodir="${genodir}"

# file paths
bedfile_pfx="${genodir}/sage" # qsub script assumes form ${bedfile_pfx}${SGE_TASK_ID}.${bedfile_sfx} 
bedfile_sfx=""
weights_file_GTEx_v6p="${dbdir}/TW_Whole_Blood_0.5_1KG.db"
weights_file_GTEx_v7="${dbdir}/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db"
weights_file_DGN="${dbdir}/DGN-WB_0.5.db"
weights_file_MESA_AFA="${dbdir}/AFA_imputed_10_peer_3_pcs_v2.db"
weights_file_MESA_AFHI="${dbdir}AFHI_imputed_10_peer_3_pcs_v2.db"
weights_file_MESA_ALL="${dbdir}ALL_imputed_10_peer_3_pcs__v2.db"
weights_file_MESA_CAU="${dbdir}/CAU_imputed_10_peer_3_pcs_v2.db"

# binaries
PYTHON=$(whereis python2.7| awk '{print $2}')
RSCRIPT=$(whereis RSCRIPT | awk '{print $2}')
PLINK=$(whereis plink | awk '{print $2}')

# scripts
BASH_run_predixcan="${thisdir}/qsub_run_predixcan.sh"
BASH_postprocess_predixcan="${thisdir}/qsub_postprocess_predixcan.sh"
R_format_predixcan="${thisdir}/format_plink_traw_into_predixcan.R"
R_impute_predixcan_genos="${thisdir}/impute_predixcan_genotypes.R" 
R_postprocess_predixcan="${thisdir}/postprocess_predixcan_results.R" 

# THIS IS VERY IMPORTANT
# there is no universal way to determine where the PrediXcan code is stored
# here it is assumed that it is cloned into a Git folder
# and that the Git folder is in the user's $HOME directory 
PREDIXCAN="${HOME}/Git/PrediXcan-master/Software/PrediXcan.py"

# SGE variables
memory_limit="2G"
scratch_memory="2G"
h_rt="0:29:59"

# output prefixes for PrediXcan
output_pfx_DGN="${dgndir}/DGN_predixcan_prediction"
output_pfx_GTEx_v6p="${gtexdir}/GTEx_v6p_predixcan_prediction"
output_pfx_GTEx_v7="${gtexdir}/GTEx_v7_predixcan_prediction"
output_pfx_MESA_AFA="${mesadir}/MESA_AFA_predixcan_prediction"
output_pfx_MESA_AFHI="${mesadir}/MESA_AFHI_predixcan_prediction"
output_pfx_MESA_ALL="${mesadir}/MESA_ALL_predixcan_prediction"
output_pfx_MESA_CAU="${mesadir}/MESA_CAU_predixcan_prediction"

# make directories if necessary
mkdir -p $dbdir
mkdir -p $genodir
mkdir -p $predixcanoutdir
mkdir -p $outdir
mkdir -p $dgndir
mkdir -p $gtexdir
mkdir -p $logdir
mkdir -p $mesadir

# job names
job_run_predixcan="sage.predixcan"
job_postprocess_predixcan="sage.predixcan.postprocess"

# ==========================================================================================
# run predixcan
# ==========================================================================================

# add binaries and scripts
qsub_variables="predixcan=$PREDIXCAN,python2=$PYTHON,Rscript=$RSCRIPT,R_format_predixcan=$R_format_predixcan,R_impute_predixcan_genos=${R_impute_predixcan_genos}"

# add directories, file paths, prefixes, suffixes
qsub_variables="${qsub_variables},outdir=$outdir,dgndir=$dgndir,gtexdir=$gtexdir,mesadir=$mesadir,predixcan_sfx=${predixcan_sfx},predixcan_genodir=${predixcan_genodir},bedfile_pfx=${bedfile_pfx},bedfile_sfx=${bedfile_sfx}"

# add PrediXcan weights
qsub_variables="${qsub_variables},weights_file_MESA_AFA=$weights_file_MESA_AFA,weights_file_MESA_AFHI=$weights_file_MESA_AFHI,weights_file_MESA_ALL=$weights_file_MESA_ALL,weights_file_MESA_CAU=$weights_file_MESA_CAU,weights_file_DGN=${weights_file_DGN},weights_file_GTEx_v6p=${weights_file_GTEx_v6p},weights_file_GTEx_v7=${weights_file_GTEx_v7}"

# output paths
qsub_variables="${qsub_variables},output_pfx_DGN=$output_pfx_DGN,output_pfx_GTEx_v6p=$output_pfx_GTEx_v6p,output_pfx_GTEx_v7=$output_pfx_GTEx_v7,output_pfx_MESA_AFA=$output_pfx_MESA_AFA,output_pfx_MESA_AFHI=$output_pfx_MESA_AFHI,output_pfx_MESA_ALL=$output_pfx_MESA_ALL,output_pfx_MESA_CAU=$output_pfx_MESA_CAU"

qsub -N ${job_run_predixcan} \
    -v ${qsub_variables} \
    -o $logdir \
    -e $logdir \
    -l mem_free=$memory_limit \
    -l h_rt=$h_rt \
    -t 1-22 \
    $BASH_run_predixcan

# ==========================================================================================
# postprocess predixcan results 
# ==========================================================================================

# binaries
qsub_variables="Rscript=$RSCRIPT"

# scripts
qsub_variables="${qsub_variables},R_postprocess_predixcan=$R_postprocess_predixcan"

# output paths
qsub_variables="${qsub_variables},output_pfx_DGN=$output_pfx_DGN,output_pfx_GTEx_v6p=$output_pfx_GTEx_v6p,output_pfx_GTEx_v7=$output_pfx_GTEx_v7,output_pfx_MESA_AFA=$output_pfx_MESA_AFA,output_pfx_MESA_AFHI=$output_pfx_MESA_AFHI,output_pfx_MESA_ALL=$output_pfx_MESA_ALL,output_pfx_MESA_CAU=$output_pfx_MESA_CAU"

qsub -N ${job_postprocess_predixcan} \
    -hold_jid ${job_run_predixcan} \
    -v ${qsub_variables} \
    -o ${logdir} \
    -e ${logdir} \
    -l mem_free=${memory_limit} \
    -l h_rt=${h_rt} \
    ${BASH_postprocess_predixcan}
