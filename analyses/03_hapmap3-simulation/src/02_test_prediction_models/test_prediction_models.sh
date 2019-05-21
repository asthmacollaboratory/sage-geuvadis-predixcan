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
#
# Other variables (nfolds_*, nthreads) are crossvalidation parameters fed to R.
# ==========================================================================================

# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -e  ## script will exit on error
set -u  ## script will exit if it sees an uninitialized variable
#set -x


# ==========================================================================================
# binaries
# ==========================================================================================
PYTHON2=$(whereis python2.7 | awk '{print $2}') ## default: system Python v2.7
RSCRIPT=$(whereis Rscript | awk '{print $2}') ## will default to system Rscript, if installed


# ==========================================================================================
# directories, filepaths
# ==========================================================================================
thisdir="$(dirname $(readlink -f $0))"
analysisdir="${thisdir}/../../analysis"
genotypedir="${analysisdir}/genotypes"
datadir="${analysisdir}/data"
outputdir="${analysisdir}/prediction_output"
output_text_dir="${outputdir}/prediction_output_text"
output_data_dir="${outputdir}/prediction_output_data"

simulation_joblist="${outputdir}/simulation_joblist.sh"
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
model_sizes=(1 5 10 20)
seeds=$(seq 2018 2117)
props=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9) ## prop = 1.0 is tested separately
nthreads=24  ## this is fed to the PYTHON script to run $NTHREADS jobs in parallel

# use these parameters for testing the script
#model_sizes=(10)
#seeds=(2018)
#props=(0.5)
#nthreads=1


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
#genenames=$(tail -n +2 $genelist | cut -f 2 | head -n 1)

# remove previous joblist if it exists
# then make empty new one
rm -f $simulation_joblist
touch $simulation_joblist

# add jobs to list
for gene in ${genenames[@]}; do

    # each gene has its own set of genotype files
    # note: some may not exist after QC!
    aa_file="${genotypedir}/AA.${gene}.chr22.raw"
    ceu_file="${genotypedir}/CEU.${gene}.chr22.raw"
    yri_file="${genotypedir}/YRI.${gene}.chr22.raw"

#    # guard against possible missing genotype files here
#    if [[ -e "${aa_file}" ]]; then
#
#        # loop over all random seeds
#        for seed in ${seeds[@]}; do
#
#            # first run cases where proportions vary
#            same_eqtls="FALSE"
#            for prop in ${props[@]}; do
#                for k in ${model_sizes[@]}; do
#
#                    # output files for nohup
#                    nohup_out="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.out"
#                    nohup_err="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.err"
#
#                    # write the actual command as echo, append to job list
#                    echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name \"${gene}\" --genotypes-pop1 \"${ceu_file}\" --genotypes-pop2 \"${yri_file}\" --genotypes-admix \"${aa_file}\" --output-directory \"${output_data_dir}\" --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop} > ${nohup_out} 2> ${nohup_err}" >> ${simulation_joblist}
#                done
#            done
#
#            # now do case where all pops share eQTLs
#            same_eqtls="TRUE"
#            prop="1.0" ## this is merely naming placeholder, it has no function when same_eqtls = TRUE
#            for k in ${model_sizes[@]}; do
#                nohup_out="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.out"
#                nohup_err="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.err"
#                echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name \"${gene}\" --genotypes-pop1 \"${ceu_file}\" --genotypes-pop2 \"${yri_file}\" --genotypes-admix \"${aa_file}\" --output-directory \"${output_data_dir}\" --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop} > ${nohup_out} 2> ${nohup_err}" >> ${simulation_joblist}
#            done
#        done
#    fi

    # now do cases where admixture proportion is varied
    # the AA file name changes by admix prop
    # thus, must loop through these separate from the previous two cases
    same_eqtls="FALSE"
    prop="0.5"
    k=10
    CEU_props=("0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1") ## <-- note that 0 is NOT 0.0, 1 is NOT 1.0 !!!
    YRI_props=("1" "0.9" "0.8" "0.7" "0.6" "0.5" "0.4" "0.3" "0.2" "0.1" "0") ## blame R for printing double as integer in analysis step 1

    # as before, assuming that #${CEU_props[@}} == #${YRI_props[@}} here...
    for i in ${!CEU_props[@]}; do

        CEU_prop=${CEU_props[$i]}
        YRI_prop=${YRI_props[$i]}
        aa_file="${genotypedir}/AA.${gene}_CEU${CEU_prop}_YRI${YRI_prop}.chr22.raw"

        if [[ -e "${aa_file}" ]]; then

        # loop over all random seeds
        for seed in ${seeds[@]}; do
                nohup_out="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.prop-${prop}.CEU-${CEU_prop}.YRI-${YRI_prop}.seed-${seed}.out"
                nohup_err="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.prop-${prop}.CEU-${CEU_prop}.YRI-${YRI_prop}.seed-${seed}.err"
                echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name \"${gene}\" --genotypes-pop1 \"${ceu_file}\" --genotypes-pop2 \"${yri_file}\" --genotypes-admix \"${aa_file}\" --output-directory \"${output_data_dir}\" --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop} --admix-proportion-pop1 ${CEU_prop} --admix-proportion-pop2 ${YRI_prop} > ${nohup_out} 2> ${nohup_err}" >> ${simulation_joblist}
        done

        fi
    done
done

## this code works for running model train-test-validate on a single computer or server
## this will schedule a Python process that will manage the threads
## think of it as a light (and brittle and fragile) SGE engine
nohup $PYTHON2 $PYTHON_multithread_commands --file ${simulation_joblist} --jobs $nthreads > ${outputdir}/nohup.multithread.out 2> ${outputdir}/nohup.multithread.err &

## need total number of jobs
#njobs=$(wc -l < ${simulation_joblist})
#
## set SGE variables
#h_rt="23:59:59"
#scratch_memory="1G"
#memory_limit="1G"
#logdir="${output_text_dir}"
#
## execute with SGE framework
#qsub -N "1kg.sim.expression" \
#     -v "command_file=${simulation_joblist}" \
#     -t 1-${njobs} \
#     -e "${logdir}" \
#     -o "${logdir}" \
#     -l mem_free="${memory_limit}" \
#     -l scratch="${scratch_memory}" \
#     -l h_rt="${h_rt}" \
#     ${BASH_qsub_template}
