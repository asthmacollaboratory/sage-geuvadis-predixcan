#!/usr/bin/env bash

# will use directory paths relative to where this script is located
# consequently, moving this script will break the script!
thisdir="$(dirname $(readlink -f $0))"
analysisdir="${thisdir}../../analysis"
crosspop_datadir="${analysisdir}/data/rnaseq/crosspop"
prediction_dir="${analysisdir}/results/elasticnet/crosspop"
outdir="${prediction_dir}/crosspop_results"
output_pfx="${outdir}/geuvadis.crosspop.results"

# make output dir if it doesn't exist already
mkdir -p ${outdir}

# file paths
RSCRIPT=$(whereis Rscript | awk '{print $2}')
R_compile_results="${thisdir}/geuvadis_compile_crosspop_prediction_results_onepop.R"
R_analyze_results="${thisdir}/geuvadis_crosspop_analyze_results.R"
popkey="${thisdir}../../datafiles/geuvadis.crosspop.sample.pop.key.txt"

# script variables
pops=("ceu" "tsi" "gbr" "yri" "fin")
altpops=("notceu" "nottsi" "notgbr" "notyri" "notfin")

# loop over the 5 training populations
# for each training pop, compile the prediction results
for i in $(seq 0 4); do
    pop="${pops[$i]}"
    altpop="${altpops[$i]}"

    rna_pop="${crosspop_datadir}/geuvadis.${pop}89.RPKM.invnorm.txt"
    rna_altpop="${crosspop_datadir}/geuvadis.${altpop}.RPKM.invnorm.txt"
    prediction_popdir="${prediction_dir}/${pop}89"
    predictions_pop="${prediction_popdir}/geuvadis_elasticnet_${pop}89_predictions.txt"
    predictions_altpop="${prediction_popdir}/geuvadis_elasticnet_${pop}89_predictinto_${altpop}.txt"

    # can parallelize this (in a fragile manner) with nohup
    # would need to separate the compilation of all pops (the part after this loop) into another script
    # for ease of execution and organization of code, run everything in serial for now
    #nohup_out="${outdir}/nohup.geuvadis.results.${pop}.out"
    #nohup_err="${outdir}/nohup.geuvadis.results.${pop}.err"

    #nohup $RSCRIPT $R_compile_results \
    $RSCRIPT $R_compile_results \
        --pop ${pop} \
        --altpop ${altpop} \
        --rna-pop ${rna_pop} \
        --rna-altpop ${rna_altpop} \
        --predictions-pop ${predictions_pop} \
        --predictions-altpop ${predictions_altpop} \
        --sample-pop-key ${popkey} \
        --output-prefix ${output_pfx} #\
        #> ${nohup_out} \
        #2> ${nohup_err} &
done

# now compile results for all pops and get summaries 
# list some variables first
crosspop_results_commongenes="${output_pfx}.*.commongenes.txt"
crosspop_results_commongenes_poscorr="${output_pfx}.*.commongenes.poscorr.txt"
allgenes="${output_pfx}.allpop.to.allpop.results.txt"
allgenes_summary="${output_pfx}.allpop.to.allpop.commongenes.summary.txt"
crosspop_results_predictions="${output_pfx}.*.predictinto.allpop.results.txt"
allgenes_poscorr="${output_pfx}.allpop.to.allpop.results.poscorr.txt"
allgenes_poscorr_summary="${output_pfx}.allpop.to.allpop.commongenes.summary.poscorr.txt"


# genes with
# -- measurements
# -- predictions in all 5 test pops
# -- predictions across all train pops
# number: 10,631
cat ${crosspop_results_commongenes} | sort | uniq -c | grep "5 E" | wc -l


# from this reduced gene list, get imputation results and put them in a new table
# first remove previous table, then add header, then add data

# construct table
rm -f ${allgenes}
touch ${allgenes}
echo -e "Gene\tTest_Pop\tCorrelation_Mean\tCorrelation_pval\tR2_Mean\tR2_pval\tTrain_Pop" > ${allgenes}
grep --no-filename -F -f <(cat ${crosspop_results_commongenes} | sort | uniq -c | grep "5 E" | awk '{ print $2 }') ${crosspop_results_predictions} >> ${allgenes} 

# summarize R2s
$RSCRIPT -e "library(data.table); library(dplyr); allgenes = fread(\"${allgenes}\"); allgenes %>% group_by(Train_Pop, Test_Pop) %>% summarize(R2_mean = mean(R2_Mean, na.rm = T), R2_StdErr = sd(R2_Mean, na.rm = T)) %>% as.data.table %>% fwrite(x = ., file = \"${allgenes_summary}\", sep = \"\\t\", na = 'NA', quote = F)"

# genes with
# -- measurements
# -- predictions in all 5 test pops
# -- predictions across all train pops 
# -- positive correlations for ALL train-test cases
# number: 142
cat ${crosspop_results_commongenes_poscorr} | sort | uniq -c | grep "5 E" | awk '{ print $2 }' | wc -l

# from this (tiny) gene list, get the imputation results
# stuff them into a new table
rm -f ${allgenes_poscorr}
touch ${allgenes_poscorr}
echo -e "Gene\tTest_Pop\tCorrelation_Mean\tCorrelation_pval\tR2_Mean\tR2_pval\tTrain_Pop" > ${allgenes_poscorr}
grep --no-filename -F -f <(cat ${crosspop_results_commongenes_poscorr} | sort | uniq -c | grep "5 E" | awk '{ print $2 }') ${crosspop_results_predictions} >> ${allgenes_poscorr} 

$RSCRIPT -e "library(data.table); library(dplyr); allgenes = fread(\"${allgenes_poscorr}\"); allgenes %>% group_by(Train_Pop, Test_Pop) %>% summarize(R2_mean = mean(R2_Mean, na.rm = T), R2_StdErr = sd(R2_Mean, na.rm = T)) %>% as.data.table %>% fwrite(x = ., file = \"${allgenes_poscorr_summary}\", sep = \"\\t\", na = 'NA', quote = F)"


# lastly, run statistical tests on the results
$RSCRIPT $R_analyze_results 
