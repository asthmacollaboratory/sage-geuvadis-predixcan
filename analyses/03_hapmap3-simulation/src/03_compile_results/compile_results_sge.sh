#!/usr/bin/env bash
# ==========================================================================================
# coded by Kevin L. Keys (2018)
# ==========================================================================================


# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -e  ## script will exit on error
set -u  ## script will exit if it sees an uninitialized variable


# ==========================================================================================
# binaries, directories, filepaths
# ==========================================================================================
#RSCRIPT=$(whereis Rscript | awk '{print $2}') ## will default to system Rscript, if installed
RSCRIPT="/usr/bin/Rscript"

thisdir="$(dirname $(readlink -f $0))"
analysisdir="${thisdir}/../../analysis"
outputdir="${analysisdir}/prediction_output"
output_text_dir="${outputdir}/prediction_output_text"
output_data_dir="${outputdir}/prediction_output_data"
logdir="${output_text_dir}/logfiles"
plotdir="${analysisdir}/plots"
resultsdir="${analysisdir}/results"
scratchdir="${HOME}/tmp"

Rdata_files="${outputdir}/Rdata_filelist.txt"
simulation_joblist="${outputdir}/simulation_joblist.sh"


# ==========================================================================================
# external scripts
# ==========================================================================================
R_compile_results="${thisdir}/compile_all_results_sge.R"
R_plot_results="${thisdir}/plot_results_sge.R"
BASH_qsub_list_files="${thisdir}/qsub_list_rdata_files.sh"
BASH_qsub_parse_results="${thisdir}/qsub_parse_results.sh"
BASH_qsub_compile_results="${thisdir}/qsub_compile_results.sh"
BASH_qsub_plot_results="${thisdir}/qsub_plot_results.sh"

# ==========================================================================================
# script variables
# ==========================================================================================
#results_file="${resultsdir}/1kg.all.results.txt"
results_file="${resultsdir}/1kg.all.results.admixvary.txt"
plot_filetype="png"
K=(1 10 20 40)

# ==========================================================================================
# executable code
# ==========================================================================================

# make directories if they don't exist
mkdir -p ${plotdir}
mkdir -p ${resultsdir}
mkdir -p ${scratchdir}
mkdir -p ${logdir}

# set SGE variables
h_rt="23:59:59"
scratch_memory="10G"
memory_limit="10G"

# specify the grid processor architecture to use
processor_architecture="lx-amd64"


#njobs=$(wc -l ${simulation_joblist} | awk '{ print $1 }') ## <-- maintains QSUB chain, but may overschedule jobs (forgiveable since extra jobs just die)
#
#### COMMENT BELOW ###
##njobs=$(wc -l ${Rdata_files} | awk '{ print $1 }') ## <-- "correct" way to schedule, but this breaks QSUB chain, must wait for previous job to finish for this to work!!
#### COMMENT ABOVE ###
#
#
## query SGE system limits, particularly the maximum permissible number of array job tasks
#max_aj_tasks=$(qconf -sconf | grep "max_aj_tasks" | awk '{ print $NF}') ## on UCSF QB3/Wynton, this is 100k
#my_max_tasks=$(echo $((max_aj_tasks / 2))) ## slightly favor fewer array jobs with more tasks vs. more array jobs with fewer tasks
#
## how many array jobs do we need to schedule?
## the +1 is necessary since BASH int division *truncates*
## if we want 3 tasks split into 2 array jobs, an int div here only yields 1 array job
## but we need *2* array jobs (1st: tasks 1,2; 2nd: task 3)
#num_array_jobs=$(echo $((njobs / my_max_tasks + 1)))
#
## will use these indexing variables to schedule tasks appropriately
#first_task=1
#last_task=${my_max_tasks}
#last_task=$(echo $(( ${last_task} > ${njobs} ? ${njobs} : ${last_task})) )
##last_task=2 ##<-- uncomment for debugging
#
## must list all previous array jobs
## this allows us to schedule everything in correct order
#sim_jobs="sim.1kg.expression.1"
#for i in $(seq 2 ${num_array_jobs}); do
#    sim_jobs="${sim_jobs},sim.1kg.expression.${i}"
#done
#
## how many data files do we have?
## will assign 1 core per data file
## must sort in order to schedule jobs in relative order
##find ${output_data_dir} -type f -name "*.Rdata" | sort > ${Rdata_files}
#qsub -N "sim.1kg.list.Rdata.files" \
#     -hold_jid "${sim_jobs}" \
#     -v "output_data_dir=${output_data_dir},Rdata_files=${Rdata_files},scratchdir=${scratchdir}" \
#     -e "${logdir}" \
#     -o "${logdir}" \
#     -l mem_free="${memory_limit}",h_rt="${h_rt}",arch="${processor_architecture}" \
#     ${BASH_qsub_list_files}
#     #-l mem_free="${memory_limit}",scratch="${scratch_memory}",h_rt="${h_rt}",arch="${processor_architecture}" \
#
#
## each iteration of this loop schedules 1 array job
#for i in $(seq 1 ${num_array_jobs}); do ##<-- COMMENT for debugging
##i=1 ##<-- uncomment for debugging, remember to comment $first_task and $last_task updates too
#
#    # compile results with SGE framework
#    # split into multiple array jobs as workaround to task limit
#    # this parses 1 data file per task, so each task is very light
#    qsub -N "sim.1kg.parse.results.${i}" \
#         -hold_jid "sim.1kg.list.Rdata.files,sim.1kg.expression.${i}" \
#         -v "RSCRIPT=${RSCRIPT},R_compile_results=${R_compile_results},scratchdir=${scratchdir},Rdata_files=${Rdata_files}" \
#         -t ${first_task}-${last_task} \
#         -e "${logdir}" \
#         -o "${logdir}" \
#         -l mem_free="${memory_limit}",h_rt="${h_rt}",arch="${processor_architecture}" \
#         ${BASH_qsub_parse_results}
#         #-l mem_free="${memory_limit}",scratch="${scratch_memory}",h_rt="${h_rt}",arch="${processor_architecture}"
#
#    first_task=$(echo $(( i*my_max_tasks + 1)) )
#    last_task=$(echo $(( (i + 1)*my_max_tasks)) )
#
#    # readjust to keep from overscheduling jobs, using ${njobs} as scheduling cap
#    first_task=$(echo $(( ${first_task} > ${njobs} ? ${njobs} : ${first_task} )) )
#    last_task=$(echo $(( ${last_task} > ${njobs} ? ${njobs} : ${last_task} )) )
#done ##<-- COMMENT for debugging
#
## update jobs on which to hold
#for i in $(seq 1 ${num_array_jobs}); do
#    sim_jobs="${sim_jobs},sim.1kg.parse.results.${i}"
#done
##sim_jobs="a"
#
## reset SGE variables
#h_rt="23:59:59"
#scratch_memory="1G"
#memory_limit="10G"
## compile all parsed results and tidy up directories
#qsub -N "sim.1kg.compile.results" \
#     -hold_jid "${sim_jobs}" \
#     -v "results_file=${results_file},scratchdir=${scratchdir}" \
#     -e "${logdir}" \
#     -o "${logdir}" \
#     -l mem_free="${memory_limit}",h_rt="${h_rt}",arch="${processor_architecture}" \
#     ${BASH_qsub_compile_results}
#     #-l mem_free="${memory_limit}",scratch="${scratch_memory}",h_rt="${h_rt}",arch="${processor_architecture}" \


# reset SGE variables
h_rt="00:10:00"
scratch_memory="1G"
memory_limit="5G"
sim_jobs="a"
sim_jobs="${sim_jobs},sim.1kg.compile.results"

# plot results within SGE framework
# run this once per # of eQTL
for k in ${K[@]}; do
    qsub -N "sim.1kg.plot.results.k${k}" \
         -hold_jid "${sim_jobs}" \
         -v "RSCRIPT=${RSCRIPT},R_plot_results=${R_plot_results},plotdir=${plotdir},plot_filetype=${plot_filetype},results_file=${results_file},k=${k}" \
         -e "${logdir}" \
         -o "${logdir}" \
         -l mem_free="${memory_limit}",h_rt="${h_rt}",arch="${processor_architecture}" \
         ${BASH_qsub_plot_results}
         #-l mem_free="${memory_limit}",scratch="${scratch_memory}",h_rt="${h_rt}",arch="${processor_architecture}" \
done
