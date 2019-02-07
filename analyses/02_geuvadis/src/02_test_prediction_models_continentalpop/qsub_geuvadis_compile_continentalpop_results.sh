#!/usr/bin/env bash

RSCRIPT=${Rscript}

R_compile_results_continentalpop="${R_compile_results_continentalpop}"

eur373_to_eur373_path="${eur373_to_eur373_path}"
eur373_to_afr_path="${eur373_to_afr_path}"
eur278_to_eur278_path="${eur278_to_eur278_path}"
eur278_to_fin_path="${eur278_to_fin_path}"
eur278_to_afr_path="${eur278_to_afr_path}"
afr_to_eur373_path="${afr_to_eur373_path}"
afr_to_afr_path="${afr_to_afr_path}"
EUR278_ids="${eur278_ids_path}"
FIN_ids="${fin_ids_path}"
output_predictions="${output_predictions}"
output_r2="${output_r2}"
poscorr="${poscorr}"
common_genes_poscorr="${common_genes_poscorr}"

$RSCRIPT $R_compile_results_continentalpop \
    --EUR373-to-EUR373 ${eur373_to_eur373_path} \
    --EUR373-to-AFR ${eur373_to_afr_path} \
    --EUR278-to-EUR278 ${eur278_to_eur278_path} \
    --EUR278-to-FIN ${eur278_to_fin_path} \
    --EUR278-to-AFR ${eur278_to_afr_path} \
    --AFR-to-EUR373 ${afr_to_eur373_path} \
    --AFR-to-AFR ${afr_to_afr_path} \
    --output-predictions ${output_predictions} \
    --output-r2 ${output_r2} \
    --EUR278-ids ${EUR278_ids} \
    --FIN-ids ${FIN_ids} \
    --poscorr ${poscorr} \
    --common-genes-poscorr ${common_genes_poscorr}

# store old command for now
#Rscript geuvadis_compile_largepop_results.R ./eur373/geuvadis_elasticnet_eur373_genelm_predvmeas_results.txt ./eur373/geuvadis_elasticnet_eur373_predictinto_yri89_genelm_predvmeas_results.txt ./eur278/geuvadis_elasticnet_eur278_genelm_predvmeas_results.txt ./eur278/geuvadis_elasticnet_eur278_predictinto_fin95_genelm_predvmeas_results.txt ./eur278/geuvadis_elasticnet_eur278_predictinto_yri89_genelm_predvmeas_results.txt ./yri89/geuvadis_elasticnet_yri89_predictinto_eur373_genelm_predvmeas_results.txt ./yri89/geuvadis_elasticnet_yri89_predictinto_eur278_genelm_predvmeas_results.txt ./yri89/geuvadis_elasticnet_yri89_predictinto_fin95_genelm_predvmeas_results.txt ./yri89/geuvadis_elasticnet_yri89_genelm_predvmeas_results.txt geuvadis.continentalpop.results.txt geuvadis.continentalpop.r2.txt geuvadis.continentalpop.r2.poscorr.txt geuvadis.continentalpop.r2.poscorr.commongenes.txt

#    --AFR-to-EUR278 ${afr_to_eur278_path} \
#    --AFR-to-FIN ${afr_to_fin_path} \
