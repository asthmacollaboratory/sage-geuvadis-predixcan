#!/usr/bin/env bash
# ==========================================================================================
# This script parses and formats SAGE expression data for use with PrediXcan. 
#
# an example call:
#
# > ./map_eqtls.sh
#
# coded by Kevin L. Keys (2017)
# ==========================================================================================

# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -e  ## script will exit on error
set -u  ## script will exit if it sees an uninitialized variable


# ==========================================================================================
# script variables
# ==========================================================================================

# script directories
thisdir="$(dirname $(readlink -f $0))"

analysisdir="${thisdir}/../../analysis"
datadir="${analysisdir}/data"
rnaseqdir="${datadir}/rnaseq"
datafiles_dir="${thisdir}/../../datafiles"
rnaseqdir="${HOME}/gala_sage/rnaseq"
alignmentdir="${rnaseqdir}/Alignments"
codedir="${rnaseqdir}/code"
resultsdir="${analysisdir}/results"
genodir="${datadir}/genotypes"
phenodir="${datadir}/phenotypes"

# make directories if they don't already exist
mkdir -p ${analysisdir}
mkdir -p ${datadir}
mkdir -p ${rnaseqdir}

# static executable paths
RSCRIPT=$(whereis RSCRIPT | awk '{print $2}')
PLINK=$(whereis plink | awk '{print $2}')

# static file paths 
vcfdir="${genodir}/vcf"
datapfx="sage"
datasfx="vcf.gz"
phenofile="${phenodir}/sage.pheno"

# script file paths
sagegenopfx="${genodir}/sage"
newsagegenopfx="${datadir}/sage.new"
mergelist="${datadir}/mergelist.txt"
pfx="sage_39_wgs_for_rnaseq"
tpfx="${datadir}/${pfx}"
traw="${tpfx}.traw"
fam="${tpfx}.fam"
genefile="${datadir}/human_ens_GRCh37_annot.extended.txt"
readsfile="${datadir}/sage_rnaseq_reads.txt"
rpmfile="${datadir}/sage_rnaseq_rpms.txt"
texprfile_normalized="${datadir}/sage_rnaseq_expressions_INTnormalized.txt"
peerfile="${datadir}/${pfx}_PEER_covariates.txt"
dups="${tpfx}.dups"
covarfile="${tpfx}.covar"
tcovarfile="${tpfx}.tcovar"
covarfile_fastqtl="${tpfx}_fastqtl.covar"
dosagefile="${tpfx}.dosages"
pcfile="${tpfx}.eigenvec"
bedfile="${tpfx}_expression.bed"
bedfile_sorted="${tpfx}_expression_sorted.bed"
bedfile_sorted_headered="${tpfx}_expression_sorted_headered.bed"
bedfile_autosomal="${tpfx}_expression_sorted_autosomal.bed"
bedfile_gz="${bedfile}.gz"
bedfile_gz_tbi="${bedfile_gz}.tbi"
rpkmfile="${datadir}/sage_rnaseq_rpkm.txt"


vcfpfx="${datadir}/joint.1484.2015_1030.ALLCHR.39.pass.x-lowComplex.dp10-gq20.SNP.id"
vcfout_uncompressed="${vcfpfx}.vcf"
vcfout="${vcfout_uncompressed}.gz"

# script names
R_normalize_reads="${thisdir}/normalize_raw_reads_to_seq_depth.R"
R_make_expr_data="${thisdir}/make_matrixeqtl_expr_file.R"
R_make_covars="${thisdir}/make_matrixeqtl_covar_file.R"
R_run_PEER="${thisdir}/run_PEER.R"

# other script variables
nthreads=4     # this applies for each chromosome
mind="0.05"    # threshold for inclusion of missingness per sample
geno="0.05"    # threshold for inclusion of missingness per genotype
min_ac=3       # minimum number for inclusion of genotype by minor allele count
hwe="0.00001"  # purge markers whose p-value for deviation from Hardy-Weinberg equilibrium falls below this number
num_peer=15    # number of PEER factors to estimate from the expression data; these are covariates in the GTEx pipeline
npcs=3         # number of PCs to compute for covariate file


# ==========================================================================================
# toggle various parts of the pipeline here 
# ==========================================================================================

# set appropriate variable to 0 to run
# any other value turns that analysis off
run_build_genotypes=1
run_build_expressions=1
run_build_covariates=1

# start time 
echo "Start time: $(date)"


# ==========================================================================================
# build the genotype file
# ==========================================================================================

if [[ "${run_build_genotypes}" -eq "0" ]]; then

    echo "building genotype file..."

    # create list of genotype VCFs for merging with bcftools
    # this assumes that they have been tabix-indexed already
    rm -f ${mergelist}
    touch ${mergelist}
    for chr in $(seq 1 22); do
        echo "${vcfdir}/chr${chr}.SAGE_hg19-1kgP3v5_panel.mixedpop_phased.eagle_39.pass.x-lowComplex.dp10-gq20.SNP.id.vcf.gz" >> ${mergelist}
    done

    # merge away! this takes awhile 
    bcftools concat \
        --file-list ${mergelist} \
        --allow-overlaps \
        --remove-duplicates \
        --output ${vcfout_uncompressed} \
        --output-type v \
        --threads 24

    # don't forget to reindex the vcf
    bgzip -c ${vcfout_uncompressed} > ${vcfout}
    tabix -f -p vcf ${vcfout}

    # must make list of duplicated SNPs
    # will feed this list into PLINK for purging duplicate SNPs
    # this approach removes *all* duplicate SNPs, taken from PLINK website:
    # https://www.cog-genomics.org/plink/1.9/data
    zcat ${vcfout} | grep -v "^#" | cut -f 3 | sort | uniq -d > ${dups}

    # with the unified VCF in hand,
    # use PLINK to convert it to BED format
    $PLINK \
        --vcf $vcfout \
        --keep-allele-order \
        --vcf-idspace-to ":" \
        --biallelic-only strict list\
        --double-id \
        --exclude ${dups} \
        --min-ac ${min_ac} \
        --hwe ${hwe} midp \
        --mind ${mind} \
        --geno ${geno} \
        --make-bed \
        --out ${tpfx} \
        --threads $nthreads

    echo -e "done building genotypes.\n"
fi

# ==========================================================================================
# build the expression file
# ==========================================================================================

if [[ "${run_build_expressions}" -eq "0" ]]; then

    echo "building expression file..."

    # if data directory does not exist then make it
    if [ ! -d $datadir ]; then
        mkdir -p $datadir;
    fi

    # remove old outfiles if present 
    rm -f $rpmfile 
    rm -f $readsfile

    # initialize new outfiles
    touch $rpmfile 
    touch $readsfile

    # get list of directories to search
    dirs=$(ls -l -d  ${alignmentdir}/*/ | awk -F " " '{ print $NF }' | sed -e 's/\/\///' | awk -F "/" '{ print $(NF-1) }')

    # need to put header to the outfile
    # peel off a column from the read file from the first directory
    genenames=$(tail -n +5 $alignmentdir/Sample1/ReadsPerGene.out.tab | awk -F " " '{ print $1 }') 
    header=$(echo $genenames)
    echo "SubjectID ${header}" > $rpmfile
    echo "SubjectID ${header}" > $readsfile

    # loop through directories
    # in each directory, pull "reads per gene" from each output file
    # clip the "//" from each directory name
    # and then append it to the outfile
    for mydir in ${dirs[@]}; do 
        fpath=$alignmentdir/$mydir/ReadsPerGene.out.tab
        reeds=$(tail -n +5 $fpath | awk -F " " '{ print $NF }') 
        rpms=$(RSCRIPT $R_normalize_reads $fpath)
        echo $mydir $reeds >> $readsfile 
        echo $mydir $rpms >> $rpmfile 
    done

    # do we need idoutfile? 
    idoutfile="${datadir}/keep_ids.txt"

    # create RPKM-normalized and transposed expression files for MatrixEQTL
    $RSCRIPT ${R_make_expr_data} \
        --RPM-file ${rpmfile} \
        --output-file ${texprfile_normalized} \
        --ID-outfile ${idoutfile} \
        --genefile ${genefile} \
        --RPKM-file ${rpkmfile} \
        --reads-file ${readsfile}

    # beat expression file into a UCSC BED
    # template taken from here:
    # https://stackoverflow.com/questions/10929453/read-a-file-line-by-line-assigning-the-value-to-a-variable
    # template reads a file line by line
    # then we extract gene name, which is 1st column
    # use the gene name to simultaneously query via grep both the gene file and the expression file
    # BED format is:
    # chr, start, end, geneID, sample1, sample2, ..., sample39
    # that is cols 5-7 of gene file, plus entire expression file
    # oddly enough, we need to write this into the header line
    # append that to ${bedfile} after writing it but before compressing it
    while read -r -a line;
    do
        # the gene name is the 1st column of the gene file
        gene="${line[0]}"

        # exit status of grep determines match
        # 0: match
        # 1: no match
        grep -q ${gene} ${texprfile_normalized}
        infile=$?

        # operate only on the matches
        if [ "${infile}" == "0" ]; then

            # this (almost) formats a line of a BED file 
            # it grabs chr, start, end from $genefile and geneID + expressions from expression file
            bedline=$(paste <(grep "${gene}" ${genefile} | cut -f 5-7 | tr "\t" " ") <( grep "${gene}" ${texprfile_normalized}));

            # change chr from just # to "chr#"
            # make sure to use TAB delimiters!
            #echo $bedline | sed -e "s/^/chr/" | tr " " "\t"
            echo $bedline | tr " " "\t"
        fi
    done < ${genefile} > ${bedfile} 

    # the BED file needs to be sorted carefully
    # meanings of arguments to sort:
    # -k1,1V: Sort 1st field alphabetically, recognising the 10 (NB -V is available in GNU coreutils >8.17) 
    # -k2,2n: Sort 2nd field numerically, loci which start first in a chromosome come first
    # -k3,3n: Sort 3rd field numerically, loci which end first come first when they have the same start position.
    cat ${bedfile} | sort -k1,1V -k2,2n -k3,3n > ${bedfile_sorted}

    # add header line
    # this line looks like
    #"#Chromosome_Name Start_Position End_Position Gene Sample1 Sample2 ..." with TAB delimiters
    # also want to purge nonautosomal expression levels
    # can do this by grepping lines which begin with a number
    headerline="#$(echo "$(head -n 1 $genefile | cut -f 5-7) $(head -n 1 $texprfile_normalized)" | tr "\t" " " | tr " " "\t")"

    # need to replace the SubjectIDs of $headerline with 
    #echo "$headerline"$'\n'"$(cat $bedfile_sorted)" | grep -P "^\d" > ${bedfile_sorted_headered}
    autosomes=$(cat $bedfile_sorted | grep -P "^\d")
    echo "${headerline}"$'\n'"${autosomes}" > ${bedfile_sorted_headered}

    echo -e "done building expression file.\n"

fi

# ==========================================================================================
# build the covariate file
# ==========================================================================================

if [[ "${run_build_covariates}" -eq "0" ]]; then

    echo "building covariate file..."

    # compute the first 3 genotype PCs as covariates in the eQTL mapping
    $PLINK \
        --bfile ${tpfx} \
        --pca ${npcs} \
        --out ${tpfx} \
        --threads ${nthreads}

    # compute PEER factors
    # GTEx uses these as covariates
    $RSCRIPT ${R_run_PEER} \
        --output_dir ${datadir} \
        --expression_file ${texprfile_normalized} \
        --prefix ${pfx} \
        --num_factors ${num_peer}

    # parse the phenotype file for covariates
    $RSCRIPT ${R_make_covars} \
        --phenotype-file ${phenofile} \
        --PC-file ${pcfile} \
        --output-file ${covarfile} \
        --PEER-file ${peerfile} \
        --transposed-output-file ${tcovarfile}

    echo -e "done.\n"

fi




# end time 
echo "Stop time: $(date)"
