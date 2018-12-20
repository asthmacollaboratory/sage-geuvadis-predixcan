#!/usr/bin/env bash
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
# This script parses GEUVADIS expression and genotype data. It produces the following files:
#
#   (1) headered and labeled matrices of expression levels, plus their transposes,
#       and lists of sample IDs, for each of EUR and YRI populations
#   (2) BGZIP'd and tabix-indexed VCFs of the gEUVADIS genotype VCFs
#   (3) a filtered, postprocessed trio of PLINK BED/BIM/FAM that covers all chromosomes
#
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

analysisdir="${thisdir}/../../analysis"
datadir="${analysisdir}/data"
geuvadis_genodir="${datadir}/genotypes"
vcfdir="${geuvadis_genodir}/vcf"
plinkdir="${geuvadis_genodir}/plink"
datafiles_dir="${thisdir}../../datafiles"
snpid_dir="${geuvadis_genodir}/snp_ids"


# ==========================================================================================
# binaries
# ==========================================================================================

# these are paths to static executables
# note that these are GUESSED with "whereis". this is a *fragile* way to set these variables!
# "whereis" can return an empty result. the guess also uses the first result, which may not be desired.
# a better solution, if it is available, is to point these variables directly to desired binaries
RSCRIPT=$(whereis Rscript | awk '{print $2}')
PLINK=$(whereis plink | awk '{print $2}')  ## <-- note: by default, doesn't use PLINK2
BCFTOOLS=$(whereis bcftools| awk '{print $2}')

# ==========================================================================================
# external scripts
# ==========================================================================================

R_remap_snpids="${thisdir}/convert_geuvadis_snpids_to_rsids.R"


# ==========================================================================================
# file paths
# ==========================================================================================

# static file paths
vcfdatapfx="GEUVADIS"
vcfdatasfx="PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes"

# script file paths
mergelist="${geuvadis_genodir}/plink/geuvadis_mergelist_for_plink.txt"
snpid_gzfile="${snpid_dir}/Phase1.Geuvadis_dbSnp137_idconvert.txt.gz"
snpid_map_file="${snpid_dir}/Phase1.Geuvadis_dbSnp137_idconvert.txt"


# ==========================================================================================
# script variables
# ==========================================================================================
nthreads=12     # this applies for each chromosome
mind="0.05"     # threshold for inclusion of missingness per sample
geno="0.05"     # threshold for inclusion of missingness per genotype
maf="0.01"      # threshold for minimum permissible minor allele frequency
hwe="0.00001"   # purge markers whose p-value test of deviation from Hardy-Weinberg equilibrium falls below this number
rsq="0.3"       # threshold for minimum permissible imputation R2


# ==========================================================================================
# executable code
# ==========================================================================================

# first make directories, if they do not already exist
mkdir -p ${analysisdir}
mkdir -p ${datadir}
mkdir -p ${geuvadis_genodir}
mkdir -p ${vcfdir}
mkdir -p ${plinkdir}
mkdir -p ${snpid_dir}


# download genotype data
for chr in $(seq 1 22); do
    current_vcf_url="https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GEUVADIS.chr${chr}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"
    wget ${current_vcf_url} -P ${vcfdir}
done


# ensure that the VCFs are all bgzip'd and tabix'd
for file in $(ls ${geuvadis_genodir}/vcf/${vcfdatapfx}.*.vcf.gz); do
    backup_vcf="${geuvadis_genodir}/vcf/backup.vcf"
    mv -f $file ${backup_vcf}
    $BCFTOOLS view --output-type z --threads $nthreads -o $file ${backup_vcf}
    ~/bin/htsfile $file;
    $BCFTOOLS index --tbi $file --threads $nthreads;
done

# download SNP ID converter
snpid_url="https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/Phase1.Geuvadis_dbSnp137_idconvert.txt.gz"
wget ${snpid_url} -P ${snpid_dir}
gzip --stdout --decompress ${snpid_gzfile} > ${snpid_map_file}

# set up a mergelist for PLINK
rm -f $mergelist
touch $mergelist

# loop over each chromosome
for chr in $(seq 1 22); do

    # file names that depend on chr
    vcfin="${vcfdir}/${vcfdatapfx}.chr${chr}.${vcfdatasfx}.vcf.gz"
    vcfout="${vcfdir}/${vcfdatapfx}.chr${chr}.${vcfdatasfx}_rsq_${rsq}.vcf.gz"
    plinkout="${plinkdir}/${vcfdatapfx}.chr${chr}.${vcfdatasfx}.rsq_${rsq}_maf_${maf}_hwe_${hwe}_geno_${geno}"

    # save current PLINK file to merge
    echo ${plinkout} >> ${mergelist}

    # now filter by imputation R2
    ${BCFTOOLS} filter --include "INFO/RSQ > ${rsq} & INFO/VT='SNP'" --output-type z --output ${vcfout} --threads ${nthreads} ${vcfin}

    # make sure to reindex with tabix protocol
    ${BCFTOOLS} index  --threads ${nthreads} --tbi ${vcfout}

    # convert VCF to PLINK format
    # apply filters along the way
    ${PLINK} \
        --vcf ${vcfout} \
        --keep-allele-order \
        --vcf-idspace-to ":" \
        --biallelic-only strict list \
        --double-id \
        --make-bed \
        --out ${plinkout} \
        --threads ${nthreads}
    ### 10 April 2018: this VCF contains both EUR and YRI, so we should not apply MAF/HWE/missingness filters yet
    ### since these data include 1KG samples, then we may not need to QC the genotypes at all
##        --hwe ${hwe} midp \
##        --maf ${maf} \
##        --geno ${geno} \

    # backup the BIM file
    rm -f ${plinkout}.bim.backup
    cp ${plinkout}.bim ${plinkout}.bim.backup

    # replace gEUVADIS SNP IDs with rsIDs
    # this will read from the *backup* BIM file and clobber the *original* BIM
    ${RSCRIPT} ${R_remap_snpids} \
        --bim-file ${plinkout}.bim.backup \
        --SNP-ID-file ${snpid_map_file} \
        --output-file ${plinkout}.bim

done

# new PLINK output prefix
plinkout="${plinkdir}/${vcfdatapfx}.ALLCHR.${vcfdatasfx}.rsq_0.8_maf_${maf}_hwe_${hwe}_geno_${geno}"

# now merge all BEDs into a single PLINK set
${PLINK} \
    --merge-list ${mergelist} \
    --make-bed \
    --out ${plinkout} \
    --threads ${nthreads}
