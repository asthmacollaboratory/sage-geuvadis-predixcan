#!/usr/bin/env bash                # -- what is the language of this shell?
#                                  # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    # -- the shell for the job                                                                               
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted                                       
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined                                         
#$ -l arch=linux-x64               # -- SGE resources (CPU type)                                                                            
# ==========================================================================================                                                
# coded by Kevin L. Keys (2019)                                                                                                             
#                                                                                                                                           
# This script prepares and executes training and testing of prediction models using                                                         
# simulated gene expression data using 1000 Genomes populations.
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

# other variables
command_file=${command_file}

# each line of ${command_file} is 1 sim configuration
# index it with (1-indexed) ${SGE_TASK_ID}, must translate to 0-index
line_number=$(expr ${SGE_TASK_ID} - 1) || true  ## guard against exit status 0

# ==========================================================================================                                                
# executable code 
# ==========================================================================================   

# start by noting current date, time, and process
echo "Date: $(date)"
echo "Host name: $(hostname)"

# which line number are we using?
echo "Running job at ${command_file}:${line_number}..."

# execute sim configuration
# note that each sim config ships with various configurations already written:
# -- correct Rscript path
# -- input paths
# -- sim parameter values
# -- output paths
# to debug, remember to look first at specific command
command $(sed -e "${$line_number}q;d" ${command_file})

# append successful exit + job info
echo "job ${JOB_ID} complete! report:"
qstat -j ${JOB_ID}
