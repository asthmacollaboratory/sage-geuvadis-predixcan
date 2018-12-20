#!/usr/bin/env bash
# ==============================================================================================================================
# copyright Asthma Collboratory (2018)
# coded by Kevin L. Keys
#
# This script computes PrediXcan weights from gEUVADIS transcriptome data.
#
# Call:
#
# ./compute_new_predixcan_weights_geuvadis.sh $ALPHA
#
# where
# -- $ALPHA = "0.0", "0.5", or "1.0", used by glmnet to determine the training algorithm
# ==============================================================================================================================

# ==============================================================================================================================
# configure BASH 
# ==============================================================================================================================
# these are used for debugging any undefined variables
#set -e
#set -u


# ==============================================================================================================================
# set script variables
# ==============================================================================================================================

# script static directories
thisdir="$(dirname $(readlink -f $0))"

MYHOME="/netapp/home/kkeys"
datadir="${thisdir}/../analysis/data"
rnaseqdir="${datadir}/rnaseq"
gctadir="${datadir}/gcta"
resultsdir="${thisdir}/../analysis/results"
glmnetdir="${resultsdir}/glmnet"
logdir="${glmnetdir}/log"
codedir="${MYHOME}/gala_sage/code"
commondir="${thisdir}../common"
imputegenodir="${MYHOME}/gala_sage/genotypes/gEUVADIS"

scratchdir="/scratch/kkeys"
outdir_eur="/scratch/kkeys/${glmmethod}/genes/eur373"
outdir_eur278="/scratch/kkeys/${glmmethod}/genes/eur278"
outdir_yri="/scratch/kkeys/${glmmethod}/genes/yri89"
outdir_tsi="/scratch/kkeys/${glmmethod}/genes/tsi"
outdir_gbr="/scratch/kkeys/${glmmethod}/genes/gbr"
outdir_fin="/scratch/kkeys/${glmmethod}/genes/fin"
outdir_ceu="/scratch/kkeys/${glmmethod}/genes/ceu"
resultsdir_eur="${glmnetdir}/${glmmethod}/eur373"
resultsdir_eur278="${glmnetdir}/${glmmethod}/eur278"
resultsdir_yri="${glmnetdir}/${glmmethod}/yri89"
resultsdir_ceu="${glmnetdir}/${glmmethod}/ceu92"
resultsdir_gbr="${glmnetdir}/${glmmethod}/gbr96"
resultsdir_fin="${glmnetdir}/${glmmethod}/fin95"
resultsdir_tsi="${glmnetdir}/${glmmethod}/tsi93"
resultsdir_crosspop="${glmnetdir}/${glmmethod}/crosspop"


resultsdir_ceu89="${resultsdir_crosspop}/ceu89"
resultsdir_tsi89="${resultsdir_crosspop}/tsi89"
resultsdir_gbr89="${resultsdir_crosspop}/gbr89"
resultsdir_fin89="${resultsdir_crosspop}/fin89"

resultssubdir_eur="${resultsdir_eur}/results"
resultssubdir_eur278="${resultsdir_eur278}/results"
resultssubdir_yri="${resultsdir_yri}/results"
resultssubdir_ceu="${resultsdir_ceu}/results"
resultssubdir_gbr="${resultsdir_gbr}/results"
resultssubdir_fin="${resultsdir_fin}/results"
resultssubdir_tsi="${resultsdir_tsi}/results"
resultssubdir_ceu89="${resultsdir_ceu89}/results"
resultssubdir_gbr89="${resultsdir_gbr89}/results"
resultssubdir_fin89="${resultsdir_fin89}/results"
resultssubdir_tsii89="${resultsdir_tsi89}/results"

tmpdir="${MYHOME}/tmp"

crosspop_dir="${rnaseqdir}/crosspop"

# make output and results directories, in case they doesn't exist
mkdir -p $outdir_eur
mkdir -p $outdir_eur278
mkdir -p $outdir_yri
mkdir -p $resultsdir_eur
mkdir -p $resultsdir_eur278
mkdir -p $resultsdir_yri
mkdir -p $resultsdir_ceu
mkdir -p $resultsdir_gbr
mkdir -p $resultsdir_fin
mkdir -p $resultsdir_tsi
mkdir -p $resultssubdir_eur
mkdir -p $resultssubdir_eur278
mkdir -p $resultssubdir_yri
mkdir -p $resultssubdir_ceu
mkdir -p $resultssubdir_gbr
mkdir -p $resultssubdir_fin
mkdir -p $resultssubdir_tsi
mkdir -p $tmpdir
mkdir -p $crosspop_dir

# filepaths 
exprfile_eur="${rnaseqdir}/geuvadis.eur373.RPKM.invnorm.txt"
exprfile_eur278="${rnaseqdir}/geuvadis.eur278.RPKM.invnorm.txt"
exprfile_yri="${rnaseqdir}/geuvadis.yri89.RPKM.invnorm.txt"
exprfile_fin="${rnaseqdir}/geuvadis.fin95.RPKM.invnorm.txt"
phenofile_eur="${rnaseqdir}/geuvadis.eur373.RPKM.invnorm.pheno"
phenofile_eur278="${rnaseqdir}/geuvadis.eur278.RPKM.invnorm.pheno"
phenofile_yri="${rnaseqdir}/geuvadis.yri89.RPKM.invnorm.pheno"
phenofile_fin="${rnaseqdir}/geuvadis.fin95.RPKM.invnorm.pheno"
genelist="${rnaseqdir}/human_ens_GRCh37_genechrpos_plusmin500kb.txt"
subjectids_eur="${rnaseqdir}/geuvadis.eur373.sampleids.txt"
subjectids_yri="${rnaseqdir}/geuvadis.yri89.sampleids.txt"
subjectids_fin="${rnaseqdir}/geuvadis.fin95.sampleids.txt"
subjectids_eur278="${rnaseqdir}/geuvadis.eur278.sampleids.txt"
predictionfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictions.txt"
predictionfile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictions.txt"
predictionfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictions.txt"
predictionfile_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89.txt"
predictionfile_eur2eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_eur373.txt"
predictionfile_eur278toyri="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_yri89.txt"
predictionfile_eur278tofin="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_fin95.txt"
predictionfile_eur278toeur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_eur278.txt"
predictionfile_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373.txt"
predictionfile_yri2yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_yri89.txt"
lambdafile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_lambdas.txt"
lambdafile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_lambdas.txt"
lambdafile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_lambdas.txt"
weightsfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_weights.txt"
weightsfile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_weights.txt"
weightsfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_weights.txt"
num_pred_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_numpred.txt"
num_pred_file_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_numpred.txt"
num_pred_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_numpred.txt"
newweightsfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_weights_noNA.txt"
newweightsfile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_weights_noNA.txt"
newweightsfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_weights_noNA.txt"
out_lm_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_lm_predvmeas_results.txt"
out_lm_file_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_lm_predvmeas_results.txt"
out_lm_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_lm_predvmeas_results.txt"
out_genelm_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_genelm_predvmeas_results.txt"
out_genelm_file_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_genelm_predvmeas_results.txt"
out_genelm_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_genelm_predvmeas_results.txt"
out_lm_file_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89_lm_predvmeas_results.txt"
out_lm_file_eur278tofin="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_fin95_lm_predvmeas_results.txt"
out_lm_file_eur278toyri="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_yri89_lm_predvmeas_results.txt"
out_lm_file_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373_lm_predvmeas_results.txt"
out_genelm_file_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89_genelm_predvmeas_results.txt"
out_genelm_file_eur278toyri="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_yri89_genelm_predvmeas_results.txt"
out_genelm_file_eur278tofin="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_fin95_genelm_predvmeas_results.txt"
out_genelm_file_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373_genelm_predvmeas_results.txt"
h2file_eur="${resultsdir_eur}/geuvadis_h2_eur373.txt"
h2file_null_eur="${resultsdir_eur}/geuvadis_h2_null_eur373.txt"
h2file_eur278="${resultsdir_eur278}/geuvadis_h2_eur278.txt"
h2file_null_eur278="${resultsdir_eur278}/geuvadis_h2_null_eur278.txt"
h2file_yri="${resultsdir_yri}/geuvadis_h2_yri89.txt"
h2file_null_yri="${resultsdir_yri}/geuvadis_h2_null_yri89.txt"

# external script locations
R_compute_new_weights="${commondir}/geuvadis_compute_new_weights.R"
R_postprocess_weights="${commondir}/geuvadis_postprocess_weights.R"
R_compute_r2="${codedir}/geuvadis_gtex_compute_r2.R"
R_subsample_eur="${codedir}/geuvadis_subsample_eur373.R" 
R_predict_new_pop="${codedir}/geuvadis_predict_in_altpop.R"
R_subsample_pop="${codedir}/geuvadis_subsample_onepop.R" 

BASH_schedule_jobs="${commondir}/geuvadis_qsub_jobs.sh"
BASH_compute_weights="${commondir}/qsub_geuvadis_compute_weights.sh"
BASH_collect_weights="${commondir}/qsub_geuvadis_collect_weights.sh"
BASH_postprocess_weights="${commondir}/qsub_geuvadis_postprocess_weights.sh"

# assumes that FIESTA is cloned into $HOME/git
# edit this if necessary
PYTHON_fiesta="${HOME}/git/albi/fiesta.py"

# binaries
PLINK=$(whereis plink | awk '{print $2}')
Rscript=$(whereis Rscript | awk '{print $2}')
GCTA=$(whereis gcta64| awk '{print $2}')
PYTHON=$(whereis python2.7| awk '{print $2}')

# script variables
maf="0.01"
hwe="0.0001"
geno="0.03"
nthreads=1
memory_limit="2G"
memory_limit_mb="2000" # manually coordinate this with $memory_limit!!!
scratch_memory="2G"
h_rt="12:00:00"
discard_ratio="0.5" # desired min ratio of samples with nonmissing LOOCV predictions, used in postprocessing, e.g. 0.5 = 50%
nsamples_eur="373"
nsamples_eur278="278"
nsamples_yri="89"
nfolds_eur="10"
nfolds_eur278="10"
nfolds_yri="89"
seed=2018

### parameters for EUR89
# set variables to population-specific parameters
# point prediction, weight, lambda, expression, subject files to EUR
# set alternative population parameters to YRI89 
h_rt="23:59:59"
pop="eur89" 
nresample="100"
eur89_dir="${rnaseqdir}/eur89"
nfolds_eur89=${nfolds_yri} ## ensure that the fold number matches the number used for training AFR 
nsamples_eur89="89"

# need headers for prediction files
predictionfile_header_eur373="Gene\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101\tHG00102\tHG00103\tHG00104\tHG00105\tHG00106\tHG00108\tHG00109\tHG00110\tHG00111\tHG00112\tHG00114\tHG00115\tHG00116\tHG00117\tHG00118\tHG00119\tHG00120\tHG00121\tHG00122\tHG00123\tHG00124\tHG00125\tHG00126\tHG00127\tHG00128\tHG00129\tHG00130\tHG00131\tHG00132\tHG00133\tHG00134\tHG00135\tHG00136\tHG00137\tHG00138\tHG00139\tHG00141\tHG00142\tHG00143\tHG00145\tHG00146\tHG00148\tHG00149\tHG00150\tHG00151\tHG00152\tHG00154\tHG00155\tHG00156\tHG00157\tHG00158\tHG00159\tHG00160\tHG00171\tHG00173\tHG00174\tHG00176\tHG00177\tHG00178\tHG00179\tHG00180\tHG00181\tHG00182\tHG00183\tHG00185\tHG00186\tHG00187\tHG00188\tHG00189\tHG00231\tHG00232\tHG00233\tHG00234\tHG00235\tHG00236\tHG00238\tHG00239\tHG00240\tHG00242\tHG00243\tHG00244\tHG00245\tHG00246\tHG00247\tHG00249\tHG00250\tHG00251\tHG00252\tHG00253\tHG00255\tHG00256\tHG00257\tHG00258\tHG00259\tHG00260\tHG00261\tHG00262\tHG00263\tHG00264\tHG00265\tHG00266\tHG00267\tHG00268\tHG00269\tHG00271\tHG00272\tHG00273\tHG00274\tHG00275\tHG00276\tHG00277\tHG00278\tHG00280\tHG00281\tHG00282\tHG00284\tHG00285\tHG00306\tHG00308\tHG00309\tHG00310\tHG00311\tHG00312\tHG00313\tHG00315\tHG00319\tHG00320\tHG00321\tHG00323\tHG00324\tHG00325\tHG00326\tHG00327\tHG00328\tHG00329\tHG00330\tHG00331\tHG00332\tHG00334\tHG00335\tHG00336\tHG00337\tHG00338\tHG00339\tHG00341\tHG00342\tHG00343\tHG00344\tHG00345\tHG00346\tHG00349\tHG00350\tHG00351\tHG00353\tHG00355\tHG00356\tHG00358\tHG00359\tHG00360\tHG00361\tHG00362\tHG00364\tHG00365\tHG00366\tHG00367\tHG00369\tHG00371\tHG00372\tHG00373\tHG00375\tHG00376\tHG00377\tHG00378\tHG00379\tHG00380\tHG00381\tHG00382\tHG00383\tHG00384\tHG01334\tHG01789\tHG01790\tHG01791\tHG02215\tNA06984\tNA06985\tNA06986\tNA06989\tNA06994\tNA07037\tNA07048\tNA07051\tNA07056\tNA07346\tNA07347\tNA07357\tNA10847\tNA10851\tNA11829\tNA11830\tNA11831\tNA11832\tNA11840\tNA11843\tNA11881\tNA11892\tNA11893\tNA11894\tNA11918\tNA11920\tNA11930\tNA11931\tNA11992\tNA11993\tNA11994\tNA11995\tNA12004\tNA12005\tNA12006\tNA12043\tNA12044\tNA12045\tNA12058\tNA12144\tNA12154\tNA12155\tNA12156\tNA12234\tNA12249\tNA12272\tNA12273\tNA12275\tNA12282\tNA12283\tNA12286\tNA12287\tNA12340\tNA12341\tNA12342\tNA12347\tNA12348\tNA12383\tNA12399\tNA12400\tNA12413\tNA12489\tNA12546\tNA12716\tNA12717\tNA12718\tNA12749\tNA12750\tNA12751\tNA12760\tNA12761\tNA12762\tNA12763\tNA12775\tNA12776\tNA12777\tNA12778\tNA12812\tNA12813\tNA12814\tNA12815\tNA12827\tNA12829\tNA12830\tNA12842\tNA12843\tNA12872\tNA12873\tNA12874\tNA12889\tNA12890\tNA20502\tNA20503\tNA20504\tNA20505\tNA20506\tNA20507\tNA20508\tNA20509\tNA20510\tNA20512\tNA20513\tNA20514\tNA20515\tNA20516\tNA20517\tNA20518\tNA20519\tNA20520\tNA20521\tNA20524\tNA20525\tNA20527\tNA20528\tNA20529\tNA20530\tNA20531\tNA20532\tNA20534\tNA20535\tNA20536\tNA20537\tNA20538\tNA20539\tNA20540\tNA20541\tNA20542\tNA20543\tNA20544\tNA20581\tNA20582\tNA20585\tNA20586\tNA20588\tNA20589\tNA20752\tNA20754\tNA20756\tNA20757\tNA20758\tNA20759\tNA20760\tNA20761\tNA20765\tNA20766\tNA20768\tNA20769\tNA20770\tNA20771\tNA20772\tNA20773\tNA20774\tNA20778\tNA20783\tNA20785\tNA20786\tNA20787\tNA20790\tNA20792\tNA20795\tNA20796\tNA20797\tNA20798\tNA20799\tNA20800\tNA20801\tNA20802\tNA20803\tNA20804\tNA20805\tNA20806\tNA20807\tNA20808\tNA20809\tNA20810\tNA20811\tNA20812\tNA20813\tNA20814\tNA20815\tNA20816\tNA20819\tNA20826\tNA20828" 
#predictionfile_header_eur373="$(head -n 1 ${exprfile_eur} | sed -e 's/,/\t/g')"

predictionfile_header_yri="Gene\tNA18486\tNA18487\tNA18488\tNA18489\tNA18498\tNA18499\tNA18502\tNA18505\tNA18508\tNA18510\tNA18511\tNA18517\tNA18519\tNA18520\tNA18858\tNA18861\tNA18867\tNA18868\tNA18870\tNA18873\tNA18907\tNA18908\tNA18909\tNA18910\tNA18912\tNA18916\tNA18917\tNA18923\tNA18933\tNA18934\tNA19092\tNA19093\tNA19095\tNA19096\tNA19098\tNA19099\tNA19102\tNA19107\tNA19108\tNA19113\tNA19114\tNA19116\tNA19117\tNA19118\tNA19119\tNA19121\tNA19129\tNA19130\tNA19131\tNA19137\tNA19138\tNA19141\tNA19143\tNA19144\tNA19146\tNA19147\tNA19149\tNA19150\tNA19152\tNA19153\tNA19159\tNA19160\tNA19171\tNA19172\tNA19175\tNA19184\tNA19185\tNA19189\tNA19190\tNA19197\tNA19198\tNA19200\tNA19201\tNA19204\tNA19206\tNA19207\tNA19209\tNA19210\tNA19213\tNA19214\tNA19222\tNA19223\tNA19225\tNA19235\tNA19236\tNA19247\tNA19248\tNA19256\tNA19257"
#predictionfile_header_yri="$(head -n 1 ${exprfile_yri89} | sed -e 's/,/\t/g')"

predictionfile_header_eur278="Gene\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101\tHG00102\tHG00103\tHG00104\tHG00105\tHG00106\tHG00108\tHG00109\tHG00110\tHG00111\tHG00112\tHG00114\tHG00115\tHG00116\tHG00117\tHG00118\tHG00119\tHG00120\tHG00121\tHG00122\tHG00123\tHG00124\tHG00125\tHG00126\tHG00127\tHG00128\tHG00129\tHG00130\tHG00131\tHG00132\tHG00133\tHG00134\tHG00135\tHG00136\tHG00137\tHG00138\tHG00139\tHG00141\tHG00142\tHG00143\tHG00145\tHG00146\tHG00148\tHG00149\tHG00150\tHG00151\tHG00152\tHG00154\tHG00155\tHG00156\tHG00157\tHG00158\tHG00159\tHG00160\tHG00231\tHG00232\tHG00233\tHG00234\tHG00235\tHG00236\tHG00238\tHG00239\tHG00240\tHG00242\tHG00243\tHG00244\tHG00245\tHG00246\tHG00247\tHG00249\tHG00250\tHG00251\tHG00252\tHG00253\tHG00255\tHG00256\tHG00257\tHG00258\tHG00259\tHG00260\tHG00261\tHG00262\tHG00263\tHG00264\tHG00265\tHG01334\tHG01789\tHG01790\tHG01791\tHG02215\tNA06984\tNA06985\tNA06986\tNA06989\tNA06994\tNA07037\tNA07048\tNA07051\tNA07056\tNA07346\tNA07347\tNA07357\tNA10847\tNA10851\tNA11829\tNA11830\tNA11831\tNA11832\tNA11840\tNA11843\tNA11881\tNA11892\tNA11893\tNA11894\tNA11918\tNA11920\tNA11930\tNA11931\tNA11992\tNA11993\tNA11994\tNA11995\tNA12004\tNA12005\tNA12006\tNA12043\tNA12044\tNA12045\tNA12058\tNA12144\tNA12154\tNA12155\tNA12156\tNA12234\tNA12249\tNA12272\tNA12273\tNA12275\tNA12282\tNA12283\tNA12286\tNA12287\tNA12340\tNA12341\tNA12342\tNA12347\tNA12348\tNA12383\tNA12399\tNA12400\tNA12413\tNA12489\tNA12546\tNA12716\tNA12717\tNA12718\tNA12749\tNA12750\tNA12751\tNA12760\tNA12761\tNA12762\tNA12763\tNA12775\tNA12776\tNA12777\tNA12778\tNA12812\tNA12813\tNA12814\tNA12815\tNA12827\tNA12829\tNA12830\tNA12842\tNA12843\tNA12872\tNA12873\tNA12874\tNA12889\tNA12890\tNA20502\tNA20503\tNA20504\tNA20505\tNA20506\tNA20507\tNA20508\tNA20509\tNA20510\tNA20512\tNA20513\tNA20514\tNA20515\tNA20516\tNA20517\tNA20518\tNA20519\tNA20520\tNA20521\tNA20524\tNA20525\tNA20527\tNA20528\tNA20529\tNA20530\tNA20531\tNA20532\tNA20534\tNA20535\tNA20536\tNA20537\tNA20538\tNA20539\tNA20540\tNA20541\tNA20542\tNA20543\tNA20544\tNA20581\tNA20582\tNA20585\tNA20586\tNA20588\tNA20589\tNA20752\tNA20754\tNA20756\tNA20757\tNA20758\tNA20759\tNA20760\tNA20761\tNA20765\tNA20766\tNA20768\tNA20769\tNA20770\tNA20771\tNA20772\tNA20773\tNA20774\tNA20778\tNA20783\tNA20785\tNA20786\tNA20787\tNA20790\tNA20792\tNA20795\tNA20796\tNA20797\tNA20798\tNA20799\tNA20800\tNA20801\tNA20802\tNA20803\tNA20804\tNA20805\tNA20806\tNA20807\tNA20808\tNA20809\tNA20810\tNA20811\tNA20812\tNA20813\tNA20814\tNA20815\tNA20816\tNA20819\tNA20826\tNA20828"
#predictionfile_header_eur278="$(head -n 1 ${exprfile_eur278} | sed -e 's/,/\t/g')"

# ==============================================================================================================================
# start script 
# ==============================================================================================================================

# -------------------- #
echo "Start time: $(date)"

# -------------------- #
# make necessary output directories
if [[ ! -d "$logdir" ]]; then mkdir -p $logdir; fi


# how many genes are in this list?
# remember that this file has a header so line count is off by 1 
nGenes=$(wc -l $genelist | cut -f 1 -d " ")
###let "nGenes=nGenes-1"

### uncomment next line for testing and debugging the pipeline
#nGenes=100
