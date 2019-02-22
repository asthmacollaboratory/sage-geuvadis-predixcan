#!/bin/bash                # -- what is the language of this shell?
#                          # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash            # -- the shell for the job
#$ -r y                    # -- tell the system that if a job crashes, it should be restarted
#$ -j y                    # -- tell the system that the STDERR and STDOUT should be joined
#$ -l arch=linux-x64       # -- SGE resources (CPU type)

# -------------------- #
# bash configuration

set -o errexit  # -- script exits after a failed command 
set -o pipefail # -- script outputs exit status of last command that threw nonzero exit status
set -o nounset  # -- script exists if it tries to use undeclared variables
#set -o xtrace   # -- uncomment for debugging; script will trace what gets executed

# user limits: -c max size of core files created
date
hostname
ulimit -c 0 

# -------------------- #
# parse script variables from SGE

# executables 
predixcan=$predixcan
python2=$python2
Rscript=$Rscript
R_format_predixcan=$R_format_predixcan
R_impute_predixcan_genos=$R_impute_predixcan_genos 

# directories
outdir=$outdir
dgndir=$dgndir
gtexdir=$gtexdir
mesadir=$mesadir

# files or file prefixes
weights_file_DGN=$weights_file_DGN
weights_file_GTEx_v6p=$weights_file_GTEx_v6p
weights_file_GTEx_v7=$weights_file_GTEx_v7
weights_file_MESA_AFA=$weights_file_MESA_AFA
weights_file_MESA_AFHI=$weights_file_MESA_AFHI
weights_file_MESA_ALL=$weights_file_MESA_ALL
weights_file_MESA_CAU=$weights_file_MESA_CAU
bedfile_pfx=$bedfile_pfx
bedfile_sfx=$bedfile_sfx
output_pfx_DGN=$output_pfx_DGN
output_pfx_GTEx_v6p=$output_pfx_GTEx_v6p
output_pfx_GTEx_v7=$output_pfx_GTEx_v7
output_pfx_MESA_AFA=$output_pfx_MESA_AFA
output_pfx_MESA_AFHI=$output_pfx_MESA_AFHI
output_pfx_MESA_ALL=$output_pfx_MESA_ALL
output_pfx_MESA_CAU=$output_pfx_MESA_CAU

# variables
chr=$SGE_TASK_ID
predixcan_sfx=${predixcan_sfx}

# -------------------- #
# begin qsub script 

# construct file path to genotypes 
bedfile="${bedfile_pfx}${chr}.${bedfile_sfx}"

# make directories that depend on current chromosome
chrdir="${outdir}/chr${chr}"
mkdir -p $chrdir

# do same for files and prefixes
predixcan_pfx="chr${chr}.${predixcan_sfx}"
predixcan_input="${predixcan_genodir}/${predixcan_pfx}"
output_pfx_DGN="${output_pfx_DGN}_chr${chr}"
output_pfx_GTEx_v6p="${output_pfx_GTEx_v6p}_chr${chr}"
output_pfx_GTEx_v7="${output_pfx_GTEx_v7}_chr${chr}"
output_pfx_MESA_AFA="${output_pfx_MESA_AFA}_chr${chr}"
output_pfx_MESA_AFHI="${output_pfx_MESA_AFHI}_chr${chr}"
output_pfx_MESA_ALL="${output_pfx_MESA_ALL}_chr${chr}"
output_pfx_MESA_CAU="${output_pfx_MESA_CAU}_chr${chr}"
samples_file="${predixcan_genodir}/chr${chr}.${bedfile_sfx}.fam"
predixcan_imputed="${chrdir}/chr${chr}.imputed.${predixcan_sfx}.gz"


# -------------------- #
# executables in this script start below

# impute missing values in PrediXcan files
# use simple Rscript to do this
# save new data.frame to file in /scratch
$Rscript $R_impute_predixcan_genos ${predixcan_input} ${predixcan_imputed}

# execute predican
# run first with GTEx v6p weights
# then run with DGN weights, then GTEx v7 weights
# lastly, run with MESA weights (African American, AfrAm+Hispanic, ALL, Caucasian)
$python2 $predixcan \
    --predict \
    --dosages ${chrdir} \
    --dosages_prefix "chr${chr}" \
    --samples ${samples_file} \
    --weights ${weights_file_GTEx_v6p} \
    --output_prefix ${output_pfx_GTEx_v6p}

$python2 $predixcan \
    --predict \
    --dosages ${chrdir} \
    --dosages_prefix "chr${chr}" \
    --samples ${samples_file} \
    --weights ${weights_file_DGN} \
    --output_prefix ${output_pfx_DGN}

$python2 $predixcan \
    --predict \
    --dosages ${chrdir} \
    --dosages_prefix "chr${chr}" \
    --samples ${samples_file} \
    --weights ${weights_file_GTEx_v7} \
    --output_prefix ${output_pfx_GTEx_v7}

$python2 $predixcan \
    --predict \
    --dosages ${chrdir} \
    --dosages_prefix "chr${chr}" \
    --samples ${samples_file} \
    --weights ${weights_file_MESA_AFA} \
    --output_prefix ${output_pfx_MESA_AFA}

$python2 $predixcan \
    --predict \
    --dosages ${chrdir} \
    --dosages_prefix "chr${chr}" \
    --samples ${samples_file} \
    --weights ${weights_file_MESA_AFHI} \
    --output_prefix ${output_pfx_MESA_AFHI}

$python2 $predixcan \
    --predict \
    --dosages ${chrdir} \
    --dosages_prefix "chr${chr}" \
    --samples ${samples_file} \
    --weights ${weights_file_MESA_ALL} \
    --output_prefix ${output_pfx_MESA_ALL}

$python2 $predixcan \
    --predict \
    --dosages ${chrdir} \
    --dosages_prefix "chr${chr}" \
    --samples ${samples_file} \
    --weights ${weights_file_MESA_CAU} \
    --output_prefix ${output_pfx_MESA_CAU}
