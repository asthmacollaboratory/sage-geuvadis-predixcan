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
# directories, filepaths 
# ==========================================================================================

thisdir="$(dirname $(readlink -f $0))"
#outputdir="${HOME}/tmp"
outputdir="${thisdir}/../../analysis/prediction_output/twas"
output_text_dir="${outputdir}/prediction_output_text"

# ensure that output directory exists
mkdir -p ${outputdir}
mkdir -p ${output_text_dir}

# make simulation job list file
# this will house all parallel jobs to execute
simulation_joblist="${outputdir}/simulation_joblist_twas.sh"
rm -f ${simulation_joblist}
touch ${simulation_joblist}

# ==========================================================================================
# script variables
# ==========================================================================================

# variables from simulation of gene expression
K=(10) ## to test all eQTL models, use K=(01 05 10 20)
proportions_eqtls=("0.0" "0.5" "0.9")
same_eqtls="FALSE"
same_effects="TRUE"

seeds=$(seq 2018 2117)
#h2_twas="0.3"
same_twas_genes="TRUE"
same_twas_effects="TRUE"
ngenes=98
ncausal_gene_numbers=(1) ## to test other causal gene models, just add here
nsamples=1000
plottype="png"
effect_sizes=("1.0" "0.5" "0.1" "0.05" "0.01" "0.005" "0.001" "0.0005" "0.0001" "0.00005" "0.00001")
nthreads=24
twas_mean="0.0" # for phenotype environmental noise
twas_sd="0.1"   # also for phenotype environmental noise

# these parameters are useful for testing the script
#seeds=$(seq 2018 2018)
#effect_sizes=("1.0")

# ==========================================================================================
# external scripts 
# ==========================================================================================
PYTHON_multithread_commands="${thisdir}/../02_test_prediction_models/multithread_commands.py"
R_simulate_twas="${thisdir}/simulate_twas.R"


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

                    # the command below is echoed to a file
                    # the job is then scheduled by the Python script run at the end
#                    $RSCRIPT $R_simulate_twas \
#                        --results-directory ${datadir} \
#                        --output-directory ${outputdir} \
#                        --output-file-prefix ${outputfilepfx} \
#                        --num-eQTLs ${k} \
#                        --proportion-shared-eQTL ${prop_shared_eqtl} \
#                        --same-eQTLs ${same_eqtls} \
#                        --same-eQTL-effects ${same_effects} \
#                        --random-seed ${seed} \
#                        --same-TWAS-genes ${same_twas_genes} \
#                        --same-TWAS-effects ${same_twas_effects} \
#                        --num-genes ${ngenes} \
#                        --num-samples ${nsamples} \
#                        --num-causal-genes ${ncausal_genes} \
#                        --plot-type ${plottype} \
#                        --effect-size ${effect_size} \
#                        --TWAS-noise-mean ${twas_mean} \
#                        --TWAS-noise-standard-deviation ${twas_sd}
#                        ##--h2-TWAS ${h2_twas} \


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

#                    $RSCRIPT $R_simulate_twas \
#                        --results-directory ${datadir} \
#                        --output-directory ${outputdir} \
#                        --output-file-prefix ${outputfilepfx} \
#                        --num-eQTLs ${k} \
#                        --proportion-shared-eQTL ${prop_shared_eqtl} \
#                        --same-eQTLs ${same_eqtls} \
#                        --same-eQTL-effects ${same_effects} \
#                        --random-seed ${seed} \
#                        --same-TWAS-genes ${same_twas_genes} \
#                        --same-TWAS-effects ${same_twas_effects} \
#                        --num-genes ${ngenes} \
#                        --num-samples ${nsamples} \
#                        --num-causal-genes ${ncausal_genes} \
#                        --plot-type ${plottype} \
#                        --effect-size ${effect_size}
#                        --TWAS-noise-mean ${twas_mean} \
#                        --TWAS-noise-standard-deviation ${twas_sd}
#                        ##--h2-TWAS ${h2_twas} \

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

# run model train-test-validate
# this will schedule a Python process that will manage the threads
# think of it as a light (and brittle and fragile) SGE engine
nohup_out="${outputdir}/nohup.twas.multithread.out"
nohup_err="${outputdir}/nohup.twas.multithread.err"
nohup $PYTHON2 $PYTHON_multithread_commands --file ${simulation_joblist} --jobs ${nthreads} > ${nohup_out} 2> ${nohup_err} &
