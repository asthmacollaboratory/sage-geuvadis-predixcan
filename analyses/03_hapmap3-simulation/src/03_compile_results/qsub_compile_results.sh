#!/usr/bin/env bash                # -- what is the language of this shell?
#                                  # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    # -- the shell for the job
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined
#$ -cwd                            # -- start the job in the current working directory
# ==========================================================================================
# coded by Kevin L. Keys (2019)
#
# This script compiles results from training and testing of prediction models using
# simulated gene expression data using 1000 Genomes populations.
# ==========================================================================================

# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -o errexit # set -e, script will exit on error
set -o nounset # set -u, script will exit if it sees an uninitialized variable
#set -o xtrace  # set -x, script will track which command is currently running

# ==========================================================================================
# script variables (passed from QSUB command)
# ==========================================================================================

# directories
scratchdir=${scratchdir}

# other variables
results_file=${results_file}

# ==========================================================================================
# executable code
# ==========================================================================================

# get absolute path of scratch directory
# this is precaution against strange behavior from GNU "find"
scratchdir_absolutepath=$(realpath ${scratchdir})

# start by noting current date, time, and process
echo "Date: $(date)"
echo "Host name: $(hostname)"

echo "This job will parse results files at ${scratchdir_absolutepath}"
echo "Results will be saved to ${results_file}"

# need header for the file
# can grab it from 1st line of any matching *_results.txt file
one_file=$(find ${scratchdir_absolutepath} -type f -name "*_results.txt" | head -n 1)
head -n 1 ${one_file} > ${results_file}

echo "Added header to results file"
echo "Compiling all results, this may take awhile..."

find $(realpath ${scratchdir_absolutepath}) -type f -name "*_results.txt" -exec grep --ignore-case --no-filename --invert-match "gene" {} +  >> ${results_file}

echo "Results parsed to file ${results_file}"
#echo "Cleaning up ${scratchdir}..."
#find ${scratchdir_absolutepath} -maxdepth 1 -type f -name "*_results.txt" -delete

# append successful exit + job info
echo "job ${JOB_ID} complete! report:"
qstat -j ${JOB_ID}
