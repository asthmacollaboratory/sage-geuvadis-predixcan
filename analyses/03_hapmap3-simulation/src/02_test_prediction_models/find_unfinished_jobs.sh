#!/usr/bin/env bash
# ==========================================================================================
# this script prints job/out/err files for jobs with missing output
# this script is extremely time-consuming and should only be used as last resort
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

simulation_joblist="${outputdir}/simulation_joblist_noRdata.sh"
simulation_outfile_list="${outputdir}/simulation_outfile_list_noRdata.sh"
simulation_errfile_list="${outputdir}/simulation_errfile_list_noRdata.sh"
genelist="${datadir}/chr22.genelist.txt"

R_simulate_crosspop_pred="${thisdir}/simulate_crosspopulation_prediction.R"

# ==========================================================================================
# code 
# ==========================================================================================

model_sizes=(10 20 40) ## will accommodate k=1 cases separately
seeds=$(seq 2018 2117)
props=("0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9") ## prop = 1.0 is tested separately
CEU_props=("0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1") ## <-- note that 0 is NOT 0.0, 1 is NOT 1.0 !!!
YRI_props=("1" "0.9" "0.8" "0.7" "0.6" "0.5" "0.4" "0.3" "0.2" "0.1" "0") ## blame R for printing double as integer in analysis step 1

same_effects="TRUE"

nfolds_external=5
nfolds_internal=10
nfolds_parallel=1

genenames=$(tail -n +2 $genelist | cut -f 2 | head -n 100)

#gene=LARGE1
#same_eqtls=FALSE
#same_effects=TRUE
#k=20
#prop="0.0"
#ceu_prop="0.2"
#yri_prop="0.8"
#seed=2018
#Rdata_name="${output_data_dir}/${gene}_simulation_prediction_admixedpop_sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}_propsharedeQTLs${prop}_CEU${ceu_prop}_YRI${yri_prop}_seed${seed}.Rdata"

# clean simulation job lists
rm -f ${simulation_joblist} ${simulation_outfile_list} ${simulation_errfile_list}
touch ${simulation_joblist} ${simulation_outfile_list} ${simulation_errfile_list}

for gene in ${genenames[@]}; do

    # each gene has its own set of genotype files
    # note: some may not exist after QC!
    aa_file="${genotypedir}/AA.${gene}.chr22.raw"
    ceu_file="${genotypedir}/CEU.${gene}.chr22.raw"
    yri_file="${genotypedir}/YRI.${gene}.chr22.raw"
    ceu_prop="0.2"
    yri_prop="0.8"

    # guard against possible missing genotype files here
    if [[ -e "${aa_file}" ]]; then

        # loop over all random seeds
        for seed in ${seeds[@]}; do

            # first run cases where proportions vary
            # do not include k=1 case here since that wastes compute resources
            same_eqtls="FALSE"
            for prop in ${props[@]}; do
                for k in ${model_sizes[@]}; do

                    Rdata_name="${output_data_dir}/${gene}_simulation_prediction_admixedpop_sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}_propsharedeQTLs${prop}_CEU${ceu_prop}_YRI${yri_prop}_seed${seed}.Rdata"
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
                        echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name ${gene} --genotypes-pop1 ${ceu_file} --genotypes-pop2 ${yri_file} --genotypes-admix ${aa_file} --output-directory ${output_data_dir} --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop}" >> ${simulation_joblist}

                    fi
                done
            done

            # k=1 case here
            k=1
            prop="0"

            Rdata_name="${output_data_dir}/${gene}_simulation_prediction_admixedpop_sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}_propsharedeQTLs${prop}_CEU${ceu_prop}_YRI${yri_prop}_seed${seed}.Rdata"
            if [[ ! -f "${Rdata_name}" ]]; then 

                nohup_out="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.out"
                nohup_err="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.err"
                echo "${nohup_out}" >> ${simulation_outfile_list}
                echo "${nohup_err}" >> ${simulation_errfile_list}
                echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name ${gene} --genotypes-pop1 ${ceu_file} --genotypes-pop2 ${yri_file} --genotypes-admix ${aa_file} --output-directory ${output_data_dir} --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop}" >> ${simulation_joblist}
            fi


            # now do case where all pops share eQTLs
            # include k=1 case here
            same_eqtls="TRUE"
            prop="1" ## this is merely naming placeholder, it has no function when same_eqtls = TRUE
            for k in 1 "${model_sizes[@]}"; do

                Rdata_name="${output_data_dir}/${gene}_simulation_prediction_admixedpop_sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}_propsharedeQTLs${prop}_CEU${ceu_prop}_YRI${yri_prop}_seed${seed}.Rdata"
                if [[ ! -f "${Rdata_name}" ]]; then 
                    nohup_out="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.out"
                    nohup_err="${output_text_dir}/${gene}.model-${same_eqtls}.effects-${same_effects}.k-${k}.seed-${seed}.prop-${prop}.err"
                    echo "$RSCRIPT $R_simulate_crosspop_pred --gene-name ${gene} --genotypes-pop1 ${ceu_file} --genotypes-pop2 ${yri_file} --genotypes-admix ${aa_file} --output-directory ${output_data_dir} --same-eqtls ${same_eqtls} --same-eqtl-effects ${same_effects} --nfolds-internal ${nfolds_internal} --nfolds-external ${nfolds_external} --nfolds-parallel ${nfolds_parallel} --num-eqtls ${k} --random-seed ${seed} --fraction-overlapping-eqtls ${prop}" >> ${simulation_joblist}
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

                Rdata_name="${output_data_dir}/${gene}_simulation_prediction_admixedpop_sameeQTLs${same_eqtls}_sameeffects${same_effects}_k${k}_propsharedeQTLs${prop}_CEU${ceu_prop}_YRI${yri_prop}_seed${seed}.Rdata"
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
