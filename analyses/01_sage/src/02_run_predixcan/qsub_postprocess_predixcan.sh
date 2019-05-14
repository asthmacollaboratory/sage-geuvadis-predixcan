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
Rscript=$Rscript
R_postprocess_predixcan=$R_postprocess_predixcan

# file path variables
output_pfx_GTEx_v6p=$output_pfx_GTEx_v6p
output_pfx_GTEx_v7=$output_pfx_GTEx_v7
output_pfx_MESA_AFA=$output_pfx_MESA_AFA
output_pfx_MESA_AFHI=$output_pfx_MESA_AFHI
output_pfx_MESA_ALL=$output_pfx_MESA_ALL
output_pfx_MESA_CAU=$output_pfx_MESA_CAU

# -------------------- #
# begin qsub script 
# executables in this script start below

# this script loops through all 22 autosomal chr files from PrediXcan
# run script one for each set of prediction weights
$Rscript $R_postprocess_predixcan $output_pfx_DGN
$Rscript $R_postprocess_predixcan $output_pfx_GTEx_v6p
$Rscript $R_postprocess_predixcan $output_pfx_GTEx_v7
$Rscript $R_postprocess_predixcan $output_pfx_MESA_AFA
$Rscript $R_postprocess_predixcan $output_pfx_MESA_AFHI
$Rscript $R_postprocess_predixcan $output_pfx_MESA_ALL
$Rscript $R_postprocess_predixcan $output_pfx_MESA_CAU
