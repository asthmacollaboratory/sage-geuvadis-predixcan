#!/usr/bin/env bash
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
# This script parses GEUVADIS gene expression data. It produces the following files:
#
#   (1) headered and labeled matrices of expression levels, plus their transposes,
#       and lists of sample IDs, for each of EUR and YRI populations
# ==========================================================================================

# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -e  ## script will exit on error
set -u  ## script will exit if it sees an uninitialized variable

# ==========================================================================================
# directories
# ==========================================================================================
thisdir="$(dirname $(readlink -f $0))"

analysisdir="${thisdir}/../analysis"
datadir="${analysisdir}/data"
geuvadis_rnaseqdir="${datadir}/rnaseq"
datafiles_dir="${thisdir}/../../datafiles"


# ==========================================================================================
# binaries
# ==========================================================================================
# these are paths to static executables
# note that these are GUESSED with "whereis". this is a *fragile* way to set these variables!
# "whereis" can return an empty result. the guess also uses the first result, which may not be desired.
# a better solution, if it is available, is to point these variables directly to desired binaries
R=$(whereis R | awk '{print $2}')
RSCRIPT=$(whereis Rscript | awk '{print $2}')

# ==========================================================================================
# external scripts
# ==========================================================================================
R_parse_geu_expr="${thisdir}/parse_geuvadis_expression_data.R"


# ==========================================================================================
# file paths, variables, URLs
# ==========================================================================================
# script file paths
all_rnaseq_data_file="${geuvadis_rnaseqdir}/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt"
eur_out_file="geuvadis.eur373.RPKM.invnorm.txt"
yri_out_file="geuvadis.yri89.RPKM.invnorm.txt"
teur_out_file="geuvadis.eur373.RPKM.invnorm.transposed.txt"
tyri_out_file="geuvadis.yri89.RPKM.invnorm.transposed.txt"
eur_ids_out_file="geuvadis.eur373.ids.txt"
yri_ids_out_file="geuvadis.yri89.ids.txt"

eur_ids_file="${datafiles_dir}/geuvadis.eur373.sampleids.txt"
yri_ids_file="${datafiles_dir}/geuvadis.yri89.sampleids.txt"

rnaseq_url="https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
rnaseq_gzfile="${geuvadis_rnaseqdir}/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"


# ==========================================================================================
# executable code
# ==========================================================================================

# first make directories, if they do not already exist
mkdir -p ${analysisdir}
mkdir -p ${datadir}
mkdir -p ${geuvadis_rnaseqdir}

## download expression data
wget ${rnaseq_url} -P ${geuvadis_rnaseqdir}
gzip --stdout --decompress ${rnaseq_gzfile} > ${all_rnaseq_data_file}

# everything required for parsing gene expression data is in one script
$RSCRIPT $R_parse_geu_expr \
    --rnaseq-data-file ${all_rnaseq_data_file} \
    --output-directory ${geuvadis_rnaseqdir} \
    --EUR-sampleIDs-file ${eur_ids_file} \
    --YRI-sampleIDs-file ${yri_ids_file} \
    --EUR-rnaseq-out-file ${eur_out_file} \
    --YRI-rnaseq-out-file ${yri_out_file} \
    --transposed-EUR-rnaseq-out-file ${teur_out_file} \
    --transposed-YRI-rnaseq-out-file ${tyri_out_file} \
    --EUR-IDs-out-file ${eur_ids_out_file} \
    --YRI-IDs-out-file ${yri_ids_out_file}
