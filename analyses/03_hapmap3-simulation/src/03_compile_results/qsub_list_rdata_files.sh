#!/usr/bin/env bash                # -- what is the language of this shell?
#                                  # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    # -- the shell for the job
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined
#$ -cwd                            # -- start the job in the current working directory
# ==========================================================================================
# coded by Kevin L. Keys (2019)
#
# This script counts the number of *.Rdata files produced in the 1000 Genomes simulation.
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
output_data_dir=${output_data_dir}

# other variables
Rdata_files=${Rdata_files}
scratchdir=${scratchdir}

# ==========================================================================================
# executable code
# ==========================================================================================

# start by noting current date, time, and process
echo "Date: $(date)"
echo "Host name: $(hostname)"

echo "This job will compile a list of Rdata files"
echo "Data directory: ${output_data_dir}"
echo "Output path: ${Rdata_files}"

# make temporary directory for "sort" command
TMPDIR="${scratchdir}/scratch_list_files"
mkdir -p ${TMPDIR}

# get list of Rdata files and save to file
#find ${output_data_dir} -maxdepth 1 -type f -name "*.Rdata" -exec sort --temporary-directory=${TMPDIR} {} \+ > ${Rdata_files}
find $(realpath ${output_data_dir}) -maxdepth 1 -type f -name "*.Rdata" -fprint ${Rdata_files}

# tidy up
rm -rf ${TMPDIR}

# append successful exit + job info
echo "job ${JOB_ID} complete! report:"
qstat -j ${JOB_ID}
