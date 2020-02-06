#!/usr/bin/env bash
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
# This script simulates two populations using HAPGEN2 and HapMap3 haplotypes from CEU + YRI.
# It first downloads HAPGEN2 and the sample HapMap3 data from the IMPUTE2 website.
# The sample haplotypes are used for forward-simulations with HAPGEN2.
# Each simulated population will contain 1000 controls and 1 case.
# The distinction is irrelevant for our purposes, and we will simply use controls.
# The simulated haplotypes from controls are eventually collapsed into genotype files.
# The script sample_CEU_YRI_haps.R constructs an admixed population from the CEU + YRI,
# and then beats everything into genotype files in PLINK RAW format.
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
download_dir="${analysisdir}/download"
bin_dir="${analysisdir}/bin"
datadir="${analysisdir}/data"
simulationdir="${analysisdir}/simulation_output"
genotypedir="${analysisdir}/genotypes"


# ==========================================================================================
# binaries
# ==========================================================================================
HAPGEN2="${bin_dir}/hapgen2"
RSCRIPT=$(whereis Rscript | awk '{print $2}') ## will default to system Rscript, if installed


# ==========================================================================================
# external scripts
# ==========================================================================================
R_sample_haps="${thisdir}/sample_CEU_YRI_haps.R"
R_get_genes="${thisdir}/get_chr22_genes.R"


# ==========================================================================================
# file paths
# ==========================================================================================
genetic_map="${datadir}/HM3/genetic_map_chr22_combined_b36.txt"
chr22_legend="${datadir}/HM3/hapmap3.r2.b36.chr22.legend"
genelist="${datadir}/chr22.genelist.txt"
CEU_hap="${datadir}/HM3/CEU.chr22.hap"
YRI_hap="${datadir}/HM3/YRI.chr22.hap"

CEU_out="${simulationdir}/CEU.chr22.out"
CEU_controls="${simulationdir}/CEU.chr22.out.controls.haps"
CEU_samples="${simulationdir}/CEU.chr22.out.controls.sample"
YRI_out="${simulationdir}/YRI.chr22.out"
YRI_controls="${simulationdir}/YRI.chr22.out.controls.haps"
YRI_samples="${simulationdir}/YRI.chr22.out.controls.sample"


# ==========================================================================================
# script variables, URLs
# ==========================================================================================
hapgen2_url="http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/download/builds/x86_64/v2.2.0/hapgen2_x86_64.tar.gz"
hm3_url="https://mathgen.stats.ox.ac.uk/wtccc-software/HM3.tgz"
hapgen2_tarball="${download_dir}/hapgen2_x86_64.tar.gz"
hm3_tarball="${download_dir}/HM3.tgz"

nCEU_haps=10000
nYRI_haps=10000


# ==========================================================================================
# executable code
# ==========================================================================================

# make directories as needed
mkdir -p ${analysisdir}
mkdir -p ${download_dir}
mkdir -p ${bin_dir}
mkdir -p ${datadir}
mkdir -p ${genotypedir}
mkdir -p ${simulationdir}

### =======================================================================================
### the next section is intentionally commented
### it shows how HAPGEN2 was used to simulate haplotypes from CEU, YRI references
### but the procedures does not yield reproducible simulations
### to reproduce analyses, instead download simulated haplotypes directly
### =======================================================================================

## download HAPGEN2 to download directory
#wget ${hapgen2_url} -P ${download_dir}
#wget ${hm3_url} -P ${download_dir}
#
## extract HAPGEN2 to bin directory
#tar -xzvf ${hapgen2_tarball} -C ${bin_dir}
#
## extract HM3 to data directory
#tar -xzvf ${hm3_tarball} -C ${datadir}
#
## grab genes from chromosome 22
## in principle, this script works for any gene and any suitable padding around the gene
## by "pad" we mean the additional bases in the cis region around the gene
## the default here is for chr22 with 500Kb pad (note: another 500Kb is added later!)
## change these defaults with optional arguments to get_chr22_genes.R, e.g.
## --cis-add 100 (adds 100 *base pairs* not Kb)
## --chr 21
#$RSCRIPT $R_get_genes --out ${genelist}
#
## a brief note on the choice of numbers for option -Ne
## (taken from http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html)
##
##     "For autosomal chromosomes, we highly recommend the values 11418 for CEPH,
##        17469 for Yoruban and 14269 for Chinese Japanese populations."
##
## presumably CEPH is the CEU? that is what is assumed here
## information for other options in HAPGEN2 is also available at the aforementioned URL
## briefly:
## -m is the genetic map (provided by IMPUTE2/HAPGEN2)
## -l is the position legend (from IMPUTE2)
## -h is the haplotype file (from HapMap3, provided by IMPUTE2)
## -o is the output path
## -dl specifies the disease SNPs and disease risks (not relevant here)
## -n is the number of *controls* and *cases* (e.g. -n 1000 1 yields 1000 controls, 1 case)
#
## simulate EUR indivs from CEU
#$HAPGEN2 \
#    -m ${genetic_map} \
#    -l ${chr22_legend} \
#    -h ${CEU_hap} \
#    -o ${CEU_out} \
#    -dl 14560203 1 1.0 2.0 \
#    -n ${nCEU_haps} 1 \
#    -Ne 11418
#
## simulate AFR indivs from YRI
#$HAPGEN2 \
#    -m ${genetic_map} \
#    -l ${chr22_legend} \
#    -h ${YRI_hap} \
#    -o ${YRI_out} \
#    -dl 14560203 1 1.0 2.0 \
#    -n ${nYRI_haps} 1 \
#    -Ne 17469

#### download required files from cloud directory
#curl -L https://ucsf.box.com/shared/static/iylsvalqrofmga8zsa8c56loj8lwznch.haps > ${CEU_controls}
#curl -L https://ucsf.box.com/shared/static/b7kpxy5t46xg6zlswqxiuoc1d4cc9wyr.sample > ${CEU_samples}
#curl -L https://ucsf.box.com/shared/static/mz0fqoldopjegbi7q216zfyfv5sy8lld.haps > ${YRI_controls}
#curl -L https://ucsf.box.com/shared/static/48qo45rqmbx84l1u1f91ek961w9mmshi.sample > ${YRI_samples}


# sample from these output files to make AA population
# then format output into genotypes in PLINK RAW format
# we only use controls, though nothing about being "control" really matters for our purposes

# start by fixing some parameters
CEU_props=("0.0" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1.0")
YRI_props=("1.0" "0.9" "0.8" "0.7" "0.6" "0.5" "0.4" "0.3" "0.2" "0.1" "0.0")
ngenes=100
gene_pad=500000 ## this is the # of base pairs in cis-region to add to gene ends
seed=2018
nCEU=10000
nYRI=10000

# assuming that #${CEU_props[@}} == #${YRI_props[@}} here...
for i in ${!CEU_props[@]}; do

    CEU_prop=${CEU_props[$i]}
    YRI_prop=${YRI_props[$i]}

    echo "CEU proportion: ${CEU_prop}"
    echo "YRI proportion: ${YRI_prop}"
    $RSCRIPT $R_sample_haps \
        --CEU-haplotype-file ${CEU_controls} \
        --YRI-haplotype-file ${YRI_controls} \
        --CEU-sample-file ${CEU_samples} \
        --YRI-sample-file ${YRI_samples} \
        --chr22-legend-path ${chr22_legend} \
        --output-dir ${genotypedir} \
        --genelist-path ${genelist} \
        --CEU-proportion ${CEU_prop} \
        --YRI-proportion ${YRI_prop} \
        --seed ${seed} \
        --num-genes ${ngenes} \
        --pad-around-genes ${gene_pad} \
        --num-CEU ${nCEU} \
        --num-YRI ${nYRI}
done
