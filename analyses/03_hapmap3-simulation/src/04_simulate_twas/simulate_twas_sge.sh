#!/usr/bin/env bash

# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -e  ## script will exit on error
set -u  ## script will exit if it sees an uninitialized variable



# ==========================================================================================
# binaries
# ==========================================================================================

PYTHON2=$(whereis python2.7 | awk '{print $2}') ## default: system Python v2.7
RSCRIPT=$(whereis Rscript | awk '{print $2}') ## will default to system Rscript, if installed


# ==========================================================================================
# directories
# ==========================================================================================

thisdir="$(dirname $(readlink -f $0))"
outputdir="${thisdir}/../../analysis/prediction_output/twas"
output_text_dir="${outputdir}/prediction_output_text"
output_results_dir="${outputdir}/prediction_output_results"
datadir="${thisdir}/../../analysis/prediction_output/prediction_output_data"
plotdir="${outputdir}/plots"

# ensure that output directory exists
mkdir -p ${outputdir}
mkdir -p ${output_text_dir}
mkdir -p ${output_results_dir}
mkdir -p ${plotdir}


# ==========================================================================================
# external scripts 
# ==========================================================================================
BASH_qsub_template="${thisdir}/../02_test_prediction_models/qsub_template.sh"
BASH_qsub_compile_results="${thisdir}/../03_compile_results/qsub_compile_results.sh"
BASH_qsub_plot_results="${thisdir}/../03_compile_results/qsub_plot_results.sh"
R_simulate_twas="${thisdir}/simulate_twas_sge.R"
R_plot_twas_results="${thisdir}/plot_twas_results_sge.R"


# ==========================================================================================
# file paths 
# ==========================================================================================

# master file for compiling results 
#results_file="${outputdir}/twas.results.txt"
results_file="${outputdir}/twas.sim.all.results.2019-06-07.txt"

# make simulation job list file
# this will house all parallel jobs to execute
simulation_joblist="${outputdir}/simulation_joblist_twas.sh"
simulation_outfile_list="${outputdir}/simulation_outfile_list_twas.sh"
simulation_errfile_list="${outputdir}/simulation_errfile_list_twas.sh"

#rm -f ${simulation_joblist}
#touch ${simulation_joblist}
#
#rm -f ${simulation_outfile_list}
#touch ${simulation_outfile_list}
#
#rm -f ${simulation_errfile_list}
#touch ${simulation_errfile_list}


# ==========================================================================================
# script variables
# ==========================================================================================

# variables from simulation of gene expression
K=(10) ## to test all eQTL models, use K=(01 05 10 20)
proportions_eqtls=("0.0" "0.5" "0.9")
same_eqtls="FALSE"
same_effects="TRUE"
plot_filetype="png"

seeds=$(seq 2018 2117)
#h2_twas="0.3"
same_twas_genes="TRUE"
same_twas_effects="TRUE"
ngenes=98
ncausal_gene_numbers=(1) ## to test other causal gene models, just add here
nsamples=1000
plottype="png"
effect_sizes=("1.0" "0.5" "0.1" "0.05" "0.025" "0.01" "0.005" "0.001" "0.0005" "0.0001" "0.00005" "0.00001")
twas_mean="0.0" # for phenotype environmental noise
twas_sd="0.1"   # also for phenotype environmental noise

# these parameters are useful for testing the script
#seeds=$(seq 2018 2018)
#effect_sizes=("1.0")


# ==========================================================================================
# executable code 
# ==========================================================================================


for k in ${K[@]}; do
    # directory of stored results depends on variables
    datadir="${HOME}/software/hapgen2/AA_sim/AA_sim_results_2018-11-30/same_eqtls_${same_eqtls}/sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}"
    for prop_shared_eqtl in ${proportions_eqtls[@]}; do
        for ncausal_genes in ${ncausal_gene_numbers[@]}; do
            for seed in ${seeds[@]}; do
                for effect_size in ${effect_sizes[@]}; do

                    outputfilepfx="temp.twas.output.k_${k}.seed_${seed}.same_eqtls_${same_eqtls}.same_eqtl_effects_${same_effects}.prop_shared_eqtl_${prop_shared_eqtl}.same_twas_genes_${same_twas_genes}.same_twas_effects_${same_twas_effects}.ncausal_genes_${ncausal_genes}.effect_size_${effect_size}"
                    twas_output_path="${outputfilepfx}.txt"
                    rm -f ${twas_output_path}
                    touch ${twas_output_path}

                    # output files for nohup
                    nohup_out="${output_text_dir}/${outputfilepfx}.out"
                    nohup_err="${output_text_dir}/${outputfilepfx}.err"

                    # write the actual command as echo, append to job list
                    echo "$RSCRIPT ${R_simulate_twas} --results-directory ${datadir} --output-directory ${outputdir} --output-file-prefix ${outputfilepfx} --num-eQTLs ${k} --proportion-shared-eQTL ${prop_shared_eqtl} --same-eQTLs ${same_eqtls} --same-eQTL-effects ${same_effects} --random-seed ${seed} --same-TWAS-genes ${same_twas_genes} --same-TWAS-effects ${same_twas_effects} --num-genes ${ngenes} --num-samples ${nsamples} --num-causal-genes ${ncausal_genes} --plot-type ${plottype} --effect-size ${effect_size} --TWAS-noise-mean ${twas_mean} --TWAS-noise-standard-deviation ${twas_sd} > ${nohup_out} 2> ${nohup_err}" >> ${simulation_joblist}
                done
            done
        done
    done
done

## repeat but for fully shared eqtls
same_eqtls="TRUE"
same_effects="TRUE"
proportions_eqtls=(1)
for k in ${K[@]}; do
    # directory of stored results depends on variables
#    datadir="${HOME}/software/hapgen2/AA_sim/AA_sim_results_2018-11-30/same_eqtls_${same_eqtls}/sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}"
    datadir="${thisdir}/../../analysis/prediction_output/prediction_output_data"
    for prop_shared_eqtl in ${proportions_eqtls[@]}; do
        for ncausal_genes in ${ncausal_gene_numbers[@]}; do
            for seed in ${seeds[@]}; do
                for effect_size in ${effect_sizes[@]}; do

                    outputfilepfx="temp.twas.output.k_${k}.seed_${seed}.same_eqtls_${same_eqtls}.same_eqtl_effects_${same_effects}.prop_shared_eqtl_${prop_shared_eqtl}.same_twas_genes_${same_twas_genes}.same_twas_effects_${same_twas_effects}.ncausal_genes_${ncausal_genes}.effect_size_${effect_size}"
                    twas_output_path="${outputfilepfx}.txt"
                    rm -f ${twas_output_path}
                    touch ${twas_output_path}

                    # output files for nohup
                    nohup_out="${output_text_dir}/${outputfilepfx}.out"
                    nohup_err="${output_text_dir}/${outputfilepfx}.err"

                    # write the actual command as echo, append to job list
                    echo "$RSCRIPT ${R_simulate_twas} --results-directory ${datadir} --output-directory ${outputdir} --output-file-prefix ${outputfilepfx} --num-eQTLs ${k} --proportion-shared-eQTL ${prop_shared_eqtl} --same-eQTLs ${same_eqtls} --same-eQTL-effects ${same_effects} --random-seed ${seed} --same-TWAS-genes ${same_twas_genes} --same-TWAS-effects ${same_twas_effects} --num-genes ${ngenes} --num-samples ${nsamples} --num-causal-genes ${ncausal_genes} --plot-type ${plottype} --effect-size ${effect_size} --TWAS-noise-mean ${twas_mean} --TWAS-noise-standard-deviation ${twas_sd} > ${nohup_out} 2> ${nohup_err}" >> ${simulation_joblist}
                done
            done
        done
    done
done

## finally, repeat for the varying admixture proportions
same_eqtls="FALSE"
same_effects="TRUE"
prop_shared_eqtl="0.5"
k=10
ncausal_genes=1
CEU_props=("0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1") ## <-- note that 0 is NOT 0.0, 1 is NOT 1.0 !!!
YRI_props=("1" "0.9" "0.8" "0.7" "0.6" "0.5" "0.4" "0.3" "0.2" "0.1" "0") ## blame R for printing double as integer in analysis step 1

# directory of stored results depends on variables
# assuming that #${CEU_props[@}} == #${YRI_props[@}} here...
for i in ${!CEU_props[@]}; do

    CEU_prop=${CEU_props[$i]}
    YRI_prop=${YRI_props[$i]}
    for seed in ${seeds[@]}; do
        for effect_size in ${effect_sizes[@]}; do

            outputfilepfx="twas.output.k_${k}.same_eqtls_${same_eqtls}.same_eqtl_effects_${same_effects}.prop_shared_eqtl_${prop_shared_eqtl}.same_twas_genes_${same_twas_genes}.same_twas_effects_${same_twas_effects}.ncausal_genes_${ncausal_genes}.seed_${seed}.CEU_${CEU_prop}.YRI_${YRI_prop}.effect_size_${effect_size}"
            twas_output_path="${outputfilepfx}.txt"
            rm -f ${twas_output_path}
            touch ${twas_output_path}

            # output files for nohup
            nohup_out="${output_text_dir}/${outputfilepfx}.out"
            nohup_err="${output_text_dir}/${outputfilepfx}.err"

            # write output log file to list
            echo "${nohup_out}" >> ${simulation_outfile_list}

            # write error log file to list
            echo "${nohup_err}" >> ${simulation_errfile_list}

            # write the actual command as echo, append to job list
            echo "$RSCRIPT ${R_simulate_twas} --results-directory ${datadir} --output-directory ${output_results_dir} --output-file-prefix ${outputfilepfx} --num-eQTLs ${k} --proportion-shared-eQTL ${prop_shared_eqtl} --same-eQTLs ${same_eqtls} --same-eQTL-effects ${same_effects} --random-seed ${seed} --same-TWAS-genes ${same_twas_genes} --same-TWAS-effects ${same_twas_effects} --num-genes ${ngenes} --num-samples ${nsamples} --num-causal-genes ${ncausal_genes} --plot-type ${plottype} --effect-size ${effect_size} --TWAS-noise-mean ${twas_mean} --TWAS-noise-standard-deviation ${twas_sd} --CEU-admixture-proportion ${CEU_prop} --YRI-admixture-proportion ${YRI_prop}" >> ${simulation_joblist}
        done
    done
done

# set SGE variables
h_rt="00:29:59"
scratch_memory="1G"
memory_limit="5G"
logdir="${output_text_dir}"

njobs=$(wc -l ${simulation_joblist} | awk '{ print $1 }')
#njobs=2

# query SGE system limits, particularly the maximum permissible number of array job tasks
max_aj_tasks=$(qconf -sconf | grep "max_aj_tasks" | awk '{ print $NF}') ## on UCSF QB3, this is 100k
my_max_tasks=$(echo $((max_aj_tasks / 2))) ## slightly favor fewer array jobs with more tasks vs. more array jobs with fewer tasks

# how many array jobs do we need to schedule?
# the +1 is necessary since BASH int division *truncates*
# if we want 3 tasks split into 2 array jobs, an int div here only yields 1 array job
# but we need *2* array jobs (1st: tasks 1,2; 2nd: task 3)
num_array_jobs=$(echo $((njobs / my_max_tasks + 1)))

# will use these indexing variables to schedule tasks appropriately
first_task=1
last_task=${my_max_tasks}
last_task=$(echo $(( ${last_task} > ${njobs} ? ${njobs} : ${last_task})) ) ## avoids overscheduling if $last_task > $njobs
#last_task=2 ##<-- uncomment for debugging

# must list *ALL* previous jobs from step 3
# this allows us to schedule everything in correct order
sim_jobs="sim.1kg.compile.results,sim.1kg.plot.results"
for i in $(seq 1 ${num_array_jobs}); do
    sim_jobs="${sim_jobs},sim.1kg.expression.${i}"
    sim_jobs="${sim_jobs},sim.1kg.parse.results.${i}"
done

# schedule TWAS simulations
# each iteration of this loop schedules 1 array job
for i in $(seq 1 ${num_array_jobs}); do
#i=1 ##<-- uncomment for debugging

    # execute with SGE framework
    qsub -N "sim.1kg.twas.${i}" \
         -v "command_file=${simulation_joblist},outfiles=${simulation_outfile_list},errfiles=${simulation_errfile_list}" \
         -t ${first_task}-${last_task} \
         -e "${logdir}" \
         -o "${logdir}" \
         -l mem_free="${memory_limit}",scratch="${scratch_memory}",h_rt="${h_rt}" \
         ${BASH_qsub_template}

    first_task=$(echo $(( i*my_max_tasks + 1)) )
    last_task=$(echo $(( (i + 1)*my_max_tasks)) )

    # readjust to keep from overscheduling jobs, using ${njobs} as scheduling cap
    first_task=$(echo $(( ${first_task} > ${njobs} ? ${njobs} : ${first_task} )) )
    last_task=$(echo $(( ${last_task} > ${njobs} ? ${njobs} : ${last_task})) )
done ##<-- comment for debugging

# update simulation jobs
for i in $(seq 1 ${num_array_jobs}); do
    sim_jobs="${sim_jobs},sim.1kg.twas.${i}"
done

# how many data files do we have?
# will assign 1 core per data file
# must sort in order to schedule jobs in relative order
#find ${output_data_dir} -type f -name "*.Rdata" | sort > ${Rdata_files}
qsub -N "sim.1kg.compile.twas.results" \
     -hold_jid "${sim_jobs}" \
     -v "scratchdir=${output_results_dir},results_file=${results_file}" \
     -e "${logdir}" \
     -o "${logdir}" \
     -l mem_free="${memory_limit}",scratch="${scratch_memory}",h_rt="${h_rt}" \
     ${BASH_qsub_compile_results}

# update simulation jobs
for i in $(seq 1 ${num_array_jobs}); do
    sim_jobs="${sim_jobs},sim.1kg.compile.twas.results"
done

# reset SGE variables
h_rt="00:29:59"
scratch_memory="1G"
memory_limit="5G"
logdir="${output_text_dir}"

## plot TWAS results
sim_jobs="a"
k=10
qsub -N "sim.1kg.plot.twas.results" \
     -hold_jid "${sim_jobs}" \
     -v "RSCRIPT=${RSCRIPT},R_plot_results=${R_plot_twas_results},results_file=${results_file},plotdir=${plotdir},plot_filetype=${plot_filetype},k=${k}" \
     -e "${logdir}" \
     -o "${logdir}" \
     -l mem_free="${memory_limit}",scratch="${scratch_memory}",h_rt="${h_rt}" \
     ${BASH_qsub_plot_results}
