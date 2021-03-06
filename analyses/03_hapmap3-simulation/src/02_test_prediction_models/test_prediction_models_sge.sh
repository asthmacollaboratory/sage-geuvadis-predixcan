#!/usr/bin/env bash
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
# This script performs train-test-validate on multiple simulated models using HAPGEN2 pops.
# Genotypes from the HAPGEN2 populations are read in PLINK RAW format.

# The models are controlled by various parameters:
# -- "same_eqtls" = TRUE or FALSE controls eQTL positions
# -- "same_effects" = TRUE or FALSE controls eQTL effect sizes
# -- "seed" controls which and how many random seeds to use
# -- "model_sizes" controls the number(s) of causal eQTLs to include in simulated models
# -- "props" controls the proportions of shared eQTLs between all populations
# -- "CEU_prop" and "YRI_prop" control the AA admixture propotions from CEU and YRI
#
# Other variables (nfolds_*) are crossvalidation parameters fed to R.
# ==========================================================================================

# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -e  ## script will exit on error
set -u  ## script will exit if it sees an uninitialized variable


# ==========================================================================================
# binaries
# ==========================================================================================
#RSCRIPT=$(whereis Rscript | awk '{print $2}') ## will default to system Rscript, if installed
RSCRIPT="/usr/bin/Rscript"


# ==========================================================================================
# directories, filepaths
# ==========================================================================================
#thisdir="$(dirname $(readlink -f $0))"
thisdir="$(readlink -f .)"
analysisdir="${thisdir}/../../analysis"
genotypedir="${analysisdir}/genotypes"
datadir="${analysisdir}/data"
outputdir="${analysisdir}/prediction_output"
output_text_dir="${outputdir}/prediction_output_text"
output_data_dir="${outputdir}/prediction_output_data"

simulation_joblist="${outputdir}/simulation_joblist.sh"
simulation_outfile_list="${outputdir}/simulation_outfile_list.sh"
simulation_errfile_list="${outputdir}/simulation_errfile_list.sh"
genelist="${datadir}/chr22.genelist.txt"


# ==========================================================================================
# external scripts
# ==========================================================================================
PYTHON_multithread_commands="${thisdir}/multithread_commands.py"
R_simulate_crosspop_pred="${thisdir}/simulate_crosspopulation_prediction.R"
BASH_qsub_template="${thisdir}/qsub_template.sh"


# ==========================================================================================
# script variables
# ==========================================================================================
same_effects="TRUE"

# these parameters control the nested crossvalidation
# in general: more $nfolds_external or $nfolds_internal --> more compututational burden
# but more $nfolds_parallel --> more parallel throughput
nfolds_external=5
nfolds_internal=10
nfolds_parallel=1

# these are the parameters used for the manuscript
#model_sizes=(1 5 10 20)
model_sizes=(10 20 40) ## will accommodate k=1 cases separately
seeds=$(seq 2018 2117)
props=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9) ## prop = 1.0 is tested separately

# use these parameters for testing the script
#model_sizes=(10)
#seeds=(2018 2019)
#props=(0.5)

# use for varying admixture proportions
#CEU_props=("0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1") ## <-- note that 0 is NOT 0.0, 1 is NOT 1.0 !!!
#YRI_props=("1" "0.9" "0.8" "0.7" "0.6" "0.5" "0.4" "0.3" "0.2" "0.1" "0") ## blame R for printing double as integer in analysis step 1
CEU_props=("0" "0.1" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1") ## <-- note that 0 is NOT 0.0, 1 is NOT 1.0 !!!
YRI_props=("1" "0.9" "0.7" "0.6" "0.5" "0.4" "0.3" "0.2" "0.1" "0") ## blame R for printing double as integer in analysis step 1
### ^^^ note absence of CEU == 0.2, YRI == 0.8; this is tested in first batch of jobs, no need to rerun it 

# specify the grid processor architecture to use
processor_architecture="lx-amd64"


# ==========================================================================================
# executable code
# ==========================================================================================

# make simulation output directory, if it doesn't exist
mkdir -p ${outputdir}
mkdir -p ${output_text_dir}
mkdir -p ${output_data_dir}

# get list of gene names to use
# use the genelist from step 1
genenames=$(tail -n +2 $genelist | cut -f 2 | head -n 100)
#genenames=$(tail -n +2 $genelist | cut -f 2 | head -n 1) ##<-- uncomment for debugging
#genenames=("ADA2" "DGCR9" "GGT1" "MRTFA")

# remove previous joblist, outfile_list, errfile_list, if they exist
# then make empty new ones
rm -f $simulation_joblist
rm -f $simulation_outfile_list
rm -f $simulation_errfile_list

touch $simulation_joblist
touch $simulation_outfile_list
touch $simulation_errfile_list

# add jobs to list
for gene in ${genenames[@]}; do

    # for first round of sims, CEU/YRI props are fixed
    CEU_prop="0.2"
    YRI_prop="0.8"

    # each gene has its own set of genotype files
    # note: some may not exist after QC!
    ceu_file="${genotypedir}/CEU.${gene}.chr22.raw"
    yri_file="${genotypedir}/YRI.${gene}.chr22.raw"
    aa_file="${genotypedir}/AA.${gene}_CEU${CEU_prop}_YRI${YRI_prop}.chr22.raw"

    # guard against possible missing genotype files here
    if [[ -e "${aa_file}" ]]; then

        # loop over all random seeds
        for seed in ${seeds[@]}; do

            # first run cases where proportions vary
            # do not include k=1 case here since that wastes compute resources
            same_eqtls="FALSE"
            for prop in ${props[@]}; do
                for k in ${model_sizes[@]}; do

                    # name of output file
                    # if this file exists, then no need to rerun the simulation for that parameter configuration
                    Rdata_name="${output_data_dir}/${gene}_simulation_prediction_admixedpop_sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}_propsharedeQTLs${prop}_CEU${CEU_prop}_YRI${YRI_prop}_seed${seed}.Rdata"

                    if [[ ! -f "${Rdata_name}" ]]; then 
                        # output log file path
                        nohup_out="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.out"

                        # error log file path
                        nohup_err="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.err"

                        # write output log file to list
                        echo "${nohup_out}" >> ${simulation_outfile_list}

                        # write error log file to list
                        echo "${nohup_err}" >> ${simulation_errfile_list}

                        # write the actual command as echo, append to job list
                        echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name ${gene} --genotypes-pop1 ${ceu_file} --genotypes-pop2 ${yri_file} --genotypes-admix ${aa_file} --output-directory ${output_data_dir} --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop} --admix-proportion-pop1 ${CEU_prop} --admix-proportion-pop2 ${YRI_prop}" >> ${simulation_joblist}
                    fi
                done
            done

            # k=1 case here
            k=1
            prop="0.0"
            Rdata_name="${output_data_dir}/${gene}_simulation_prediction_admixedpop_sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}_propsharedeQTLs${prop}_CEU${CEU_prop}_YRI${YRI_prop}_seed${seed}.Rdata"
            if [[ ! -f "${Rdata_name}" ]]; then 
            
                nohup_out="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.out"
                nohup_err="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.err"
                echo "${nohup_out}" >> ${simulation_outfile_list}
                echo "${nohup_err}" >> ${simulation_errfile_list}
                echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name ${gene} --genotypes-pop1 ${ceu_file} --genotypes-pop2 ${yri_file} --genotypes-admix ${aa_file} --output-directory ${output_data_dir} --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop} --admix-proportion-pop1 ${CEU_prop} --admix-proportion-pop2 ${YRI_prop}" >> ${simulation_joblist}

            fi

            # now do case where all pops share eQTLs
            # include k=1 case here
            same_eqtls="TRUE"
            prop="1.0" ## this is merely naming placeholder, it has no function when same_eqtls = TRUE
            for k in 1 "${model_sizes[@]}"; do

                Rdata_name="${output_data_dir}/${gene}_simulation_prediction_admixedpop_sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}_propsharedeQTLs${prop}_CEU${CEU_prop}_YRI${YRI_prop}_seed${seed}.Rdata"
                if [[ ! -f "${Rdata_name}" ]]; then 

                    nohup_out="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.out"
                    nohup_err="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.err"
                    echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name ${gene} --genotypes-pop1 ${ceu_file} --genotypes-pop2 ${yri_file} --genotypes-admix ${aa_file} --output-directory ${output_data_dir} --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop} --admix-proportion-pop1 ${CEU_prop} --admix-proportion-pop2 ${YRI_prop}" >> ${simulation_joblist}
                    echo "${nohup_out}" >> ${simulation_outfile_list}
                    echo "${nohup_err}" >> ${simulation_errfile_list}

                fi
            done
        done
    fi

    # now do cases where admixture proportion is varied
    # the AA file name changes by admix prop
    # thus, must loop through these separate from the previous two cases
    same_eqtls="FALSE"
    prop="0.5"
    k=10

    # assuming that #${CEU_props[@}} == #${YRI_props[@}} here...
    for i in ${!CEU_props[@]}; do

        CEU_prop=${CEU_props[$i]}
        YRI_prop=${YRI_props[$i]}
        aa_file="${genotypedir}/AA.${gene}_CEU${CEU_prop}_YRI${YRI_prop}.chr22.raw"

        # guard against possible missing genotype files here
        if [[ -e "${aa_file}" ]]; then

            # loop over all random seeds
            for seed in ${seeds[@]}; do

                Rdata_name="${output_data_dir}/${gene}_simulation_prediction_admixedpop_sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}_propsharedeQTLs${prop}_CEU${CEU_prop}_YRI${YRI_prop}_seed${seed}.Rdata"
                if [[ ! -f "${Rdata_name}" ]]; then 

                    # output log file path
                    nohup_out="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.prop-${prop}.CEU-${CEU_prop}.YRI-${YRI_prop}.seed-${seed}.out"

                    # error log file path
                    nohup_err="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.prop-${prop}.CEU-${CEU_prop}.YRI-${YRI_prop}.seed-${seed}.err"

                    # write Rscript command to list
                    echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name ${gene} --genotypes-pop1 ${ceu_file} --genotypes-pop2 ${yri_file} --genotypes-admix ${aa_file} --output-directory ${output_data_dir} --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop} --admix-proportion-pop1 ${CEU_prop} --admix-proportion-pop2 ${YRI_prop}" >> ${simulation_joblist}

                    # write output log file to list
                    echo "${nohup_out}" >> ${simulation_outfile_list}

                    # write error log file to list
                    echo "${nohup_err}" >> ${simulation_errfile_list}
                fi
            done
        fi
    done
done


# set SGE variables
h_rt="00:29:59"
scratch_memory="1G"
memory_limit="10G"
logdir="${output_text_dir}/logfiles"
mkdir -p ${logdir}

# need total number of jobs
njobs=$(wc -l < ${simulation_joblist})
#njobs=20 ##<-- uncomment for debugging

# query SGE system limits, particularly the maximum permissible number of array job tasks
max_aj_tasks=$(qconf -sconf | grep "max_aj_tasks" | awk '{ print $NF}') ## on UCSF QB3/Wynton, this is 100k
#my_max_tasks=$(echo $((max_aj_tasks / 2))) ## paranoid safety precaution, keeps array jobs within limits
my_max_tasks=$(echo $((max_aj_tasks / 10))) ## favor fewer array jobs with more tasks vs. more array jobs with fewer tasks 

# how many array jobs do we need to schedule?
# the +1 is necessary since BASH int division does *truncation*
# therefore, for 3 tasks in groups of 2, the int div yields 1 array job
# but we need *2* array jobs (1st: tasks 1,2; 2nd: task 3)
num_array_jobs=$(echo $((njobs / my_max_tasks + 1)))

# will use these indexing variables to schedule tasks appropriately
first_task=1
last_task=${my_max_tasks}

# avoid overscheduling jobs if $njobs is smaller than $my_max_tasks
last_task=$(echo $(( ${last_task} > ${njobs} ? ${njobs} : ${last_task})) )
#last_task=20 ##<-- uncomment for debugging

# each iteration of this loop schedules 1 array job
for i in $(seq 1 ${num_array_jobs}); do ##<-- COMMENT HERE for debugging
#i=1 ##<-- uncomment for debugging

    array_job_name="sim.1kg.expression.${i}"
    array_job_logdir="${logdir}/${array_job_name}"
    mkdir -p ${array_job_logdir}
    # execute with SGE framework
    qsub -N "${array_job_name}" \
         -v "command_file=${simulation_joblist},outfiles=${simulation_outfile_list},errfiles=${simulation_errfile_list}" \
         -t ${first_task}-${last_task} \
         -e "${array_job_logdir}" \
         -o "${array_job_logdir}" \
         -l mem_free="${memory_limit}",h_rt="${h_rt}",arch="${processor_architecture}" \
         ${BASH_qsub_template}
         #-l scratch="${scratch_memory}" \

    first_task=$(echo $(( i*my_max_tasks + 1)) )
    last_task=$(echo $(( (i + 1)*my_max_tasks)) )

    # readjust to keep from overscheduling jobs, using ${njobs} as scheduling cap
    first_task=$(echo $(( ${first_task} > ${njobs} ? ${njobs} : ${first_task} )) )
    last_task=$(echo $(( ${last_task} > ${njobs} ? ${njobs} : ${last_task})) )
done ##<-- COMMENT HERE for debugging
