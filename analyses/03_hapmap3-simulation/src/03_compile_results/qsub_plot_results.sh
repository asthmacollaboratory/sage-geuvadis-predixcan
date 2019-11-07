#!/usr/bin/env bash                # -- what is the language of this shell?
#                                  # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    # -- the shell for the job
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined
#$ -cwd                            # -- start the job in the current working directory
#$ -l arch=linux-x64               # -- SGE resources (CPU type)
# ==========================================================================================
# coded by Kevin L. Keys (2019)
#
# This script plots results from analyses using simulated gene expression data
# using 1000 Genomes populations.
# ==========================================================================================

# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -o errexit # set -e, script will exit on error
set -o nounset # set -u, script will exit if it sees an uninitialized variable
#set -o xtrace  # set -x, for debugging, track which command is currently running

# ==========================================================================================
# script variables (passed from QSUB command)
# ==========================================================================================

# binaries
RSCRIPT=${RSCRIPT}

# scripts
R_plot_results=${R_plot_results}

# directories
plotdir=${plotdir}

# other variables
results_file=${results_file}
plot_filetype=${plot_filetype}
k=${k}

# ==========================================================================================
# executable code
# ==========================================================================================

# start by noting current date, time, and process
echo "Date: $(date)"
echo "Host name: $(hostname)"

# plot results and save to $plotdir
$RSCRIPT $R_plot_results \
    --results-file ${results_file} \
    --output-directory ${plotdir} \
    --plot-filetype ${plot_filetype} \
    --num-eQTL ${k}

# append successful exit + job info
echo "job ${JOB_ID} complete! report:"
qstat -j ${JOB_ID}
