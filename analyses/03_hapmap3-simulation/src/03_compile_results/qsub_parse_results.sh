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
# This script compiles results from training and testing of prediction models using
# simulated gene expression data using 1000 Genomes populations.
# ==========================================================================================

# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -o errexit # set -e, script will exit on error
set -o nounset # set -u, script will exit if it sees an uninitialized variable
#set -o xtrace  # set -x, script will track which command is currently running
#ulimit -c 0 # user limits: -c covers the max size of core files created

# ==========================================================================================
# script variables (passed from QSUB command)
# ==========================================================================================

# binaries
RSCRIPT=${RSCRIPT}

# scripts
R_compile_results=${R_compile_results}

# directories
scratchdir=${scratchdir}

# other variables
Rdata_files=${Rdata_files}

# each line of ${command_file} is 1 sim configuration
# index it with (1-indexed) ${SGE_TASK_ID}
line_number=${SGE_TASK_ID}

# ==========================================================================================
# executable code
# ==========================================================================================

# start by noting current date, time, and process
echo "Date: $(date)"
echo "Host name: $(hostname)"

# which line number are we using?
echo "This job will parse file at ${Rdata_files}:${line_number}..."

# set paths for data file, results file
data_file=$(sed -e "${line_number}q;d" ${Rdata_files})
results_file="$(echo ${data_file} | awk -F "/" '{ print $NF}' | sed -e 's/\.Rdata//')_results.txt"

echo "Processing file ${data_file}"

# compile results from current Rdata file and save to $results_file
$RSCRIPT $R_compile_results \
    --results-file ${data_file} \
    --output-directory ${scratchdir} \
    --output-filename ${results_file}

echo "Results parsed to file ${results_file}"

# append successful exit + job info
echo "job ${JOB_ID} complete! report:"
qstat -j ${JOB_ID}
