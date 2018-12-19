###!/usr/bin/env Rscript --vanilla
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
# This script simulates phenotypes and computes predictive models in three populations.
# The populations are specified by genotype files from which phenotypes are simulated. 
# Each phenotype contains a prespecified number of true eQTLs; other genotypes have 0 effect.
# Each phenotype/genotype set is then used to train a predictive model.
# The training scheme uses elastic net regression with a nested crossvalidation scheme.
# ==========================================================================================

# ==========================================================================================
# script options 
# ==========================================================================================

# don't emit warnings, only errors
#options(warn = -1)

# ==========================================================================================
# libraries
# ==========================================================================================
suppressMessages(library(data.table))
suppressMessages(library(methods))
suppressMessages(library(optparse))


# parse command line variables
option_list = list(
    make_option(
        c("-ch", "--CEU-haplotype-file"),
        type    = "character",
        default = NULL, 
        help    = "HAPGEN2 haplotype file for HapMap3 CEU",
        metavar = "character"
    ),
    make_option(
        c("-yh", "--YRI-haplotype-file"),
        type    = "character",
        default = NULL, 
        help    = "HAPGEN2 haplotype file for HapMap3 YRI", 
        metavar = "character"
    ),
    make_option(
        c("-cs", "--CEU-sample-file"),
        type    = "character",
        default = NULL, 
        help    = "HAPGEN2 sample file for HapMap3 CEU",
        metavar = "character"
    ),
    make_option(
        c("-ys", "--YRI-sample-file"),
        type    = "character",
        default = NULL, 
        help    = "HAPGEN2 sample file for HapMap3 YRI", 
        metavar = "character"
    ),
    make_option(
        c("-l", "--chr22-legend-path"),
        type    = "character",
        default = NULL, 
        help    = "Markers on chromosome 22, from HAPGEN2, e.g. hapmap3.r2.b36.chr22.legend", 
        metavar = "character"
    ),
    make_option(
        c("-g", "--genelist-path"),
        type    = "character",
        default = NULL, 
        help    = "List of genes on chr22 to use for simulation, e.g. 'chr22.genelist.txt'", 
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-dir"),
        type    = "character",
        default = "../../output", 
        help    = "Directory for output [default = %default]", 
        metavar = "character"
    ),
    make_option(
        c("-cp", "--CEU-proportion"),
        type    = "double",
        default = 0.2, 
        help    = "Proportion of CEU haplotypes to draw for admixed population [default = %default]",
        metavar = "character"
    ),
    make_option(
        c("-yp", "--YRI-proportion"),
        type    = "double",
        default = 0.8, 
        help    = "Proportion of YRI haplotypes to draw for admixed population [default = %default]",
        metavar = "character"
    ),
    make_option(
        c("-s", "--seed"),
        type    = "integer",
        default = 2018, 
        help    = "Random seed to fix sampling of haplotypes [default = %default]", 
        metavar = "integer"
    ),
    make_option(
        c("-p", "--pad-around-genes"),
        type    = "integer",
        default = 500000, 
        help    = "Number of base pairs to add to ends of genes [default = %default]",
        metavar = "integer"
    ),
    make_option(
        c("-n", "--num-genes"),
        type    = "integer",
        default = 100, 
        help    = "Number of genes to parse, starting from largest [default = %default]",
        metavar = "integer"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

# ==========================================================================================
# script variables
# ==========================================================================================

ceu.hap.file.path = opt$CEU_haplotype_file 
yri.hap.file.path = opt$YRI_haplotype_file 
ceu.samples.path  = opt$CEU_sample_file 
yri.samples.path  = opt$YRI_sample_file
chr22.legend.path = opt$chr22_legend_path
output.dir        = opt$output_dir
seed              = opt$seed
ceu.prop          = opt$CEU_proportion 
yri.prop          = opt$YRI_proportion 
genelist.path     = opt$genelist_path
pad               = opt$pad_around_genes
ngenes            = opt$num_genes

# set random seed
set.seed(seed)

# load haps
ceu.haps = fread(ceu.hap.file.path)
yri.haps = fread(yri.hap.file.path)

# sample proportions
ceu.num = dim(ceu.haps)[2]
yri.num = dim(yri.haps)[2]
ceu.haps.sub = ceu.haps[,sample(1:ceu.num,floor(ceu.prop*ceu.num)), with = FALSE]
yri.haps.sub = yri.haps[,sample(1:yri.num,ceiling(yri.prop*yri.num)), with = FALSE]

# slap haps together
aa.haps = data.table(cbind(ceu.haps.sub, yri.haps.sub))
aa.num  = dim(aa.haps)[2]

# write to file
aa.haps.file = file.path(output.dir, "AA.chr22.hap")
fwrite(x = aa.haps, file = aa.haps.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

# create genotype files
# do this by collapsing every other column 
ceu.geno = data.table(ceu.haps[,seq(1,ceu.num,2), with = FALSE] + ceu.haps[,seq(2,ceu.num,2), with = FALSE])
yri.geno = data.table(yri.haps[,seq(1,yri.num,2), with = FALSE] + yri.haps[,seq(2,yri.num,2), with = FALSE])
aa.geno  = data.table(aa.haps[,seq(1,aa.num,2), with = FALSE]   + aa.haps[,seq(2,aa.num,2), with = FALSE])

# the $*.geno files contain SNPs in rows and genotypes in columns
# transpose to SNP-major format (SNPs on columns
ceu.geno = data.table(t(ceu.geno))
yri.geno = data.table(t(yri.geno))
aa.geno  = data.table(t(aa.geno))

# write genotype files to file
aa.geno.file  = file.path(output.dir, "AA.chr22.geno")
ceu.geno.file = file.path(output.dir, "CEU.chr22.geno")
yri.geno.file = file.path(output.dir, "YRI.chr22.geno")

fwrite(x = aa.geno,  file = aa.geno.file,  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
fwrite(x = ceu.geno, file = ceu.geno.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
fwrite(x = yri.geno, file = yri.geno.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# slap sample information onto left of subset genotype files
# want to make a PLINK RAW file
# leftmost 6 cols: FID, IID, PAT, MAT, SEX, PHENO
# first 4 come from map file, latter two set to 0 or -9 ("missing" in PLINK)
# note: CEU samples for AA works because the samples file is mostly a placeholder
aa.samples  = fread(ceu.samples.path, header = TRUE)
ceu.samples = fread(ceu.samples.path, header = TRUE)
yri.samples = fread(yri.samples.path, header = TRUE)

aa.samples  = aa.samples[-1,]
ceu.samples = ceu.samples[-1,]
yri.samples = yri.samples[-1,]

colnames(aa.samples)  = c("FID", "IID", "PAT", "MAT")
colnames(ceu.samples) = c("FID", "IID", "PAT", "MAT")
colnames(yri.samples) = c("FID", "IID", "PAT", "MAT")

# this column is a placeholder and not needed for our purposes
aa.samples$SEX  = 0
ceu.samples$SEX = 0
yri.samples$SEX = 0

# this column is also a placeholder and not needed for our purposes
aa.samples$PHENO  = -9
ceu.samples$PHENO = -9
yri.samples$PHENO = -9

# want to parse genes from these data
# use the chr22 gene list created earlier
chr22.leg = fread(chr22.legend.path)
genelist  = fread(genelist.path)

# ensure that the requested number of genes is not more than what is available in genelist
ngenes.inlist = dim(genelist)[1]
if ( ngenes.inlist < ngenes ){
    warning(paste0("Genelist contains ", ngenes.inlist, ", but --num-genes=", ngenes, " requested. Will subset ", ngenes.inlist, " instead.\n"))
    ngenes = ngenes.inlist
}

# loop over requested number of genes in $genelist
for (i in 1:ngenes) {

    # get HGNC name of the current gene
    genename  = genelist$hgnc_symbol[i]

    # grab Boolean vector of SNPs that lie in current gene's boundary
    # argument $pad assumes that we haven't added any additional cis markers to the gene
    # note that by default $pad = 500kb, so if $genelist already padded the genes then this code by default adds more padding! 
    gene.snps = genelist$start_position[i] - pad < chr22.leg$pos & genelist$end_position[i] + pad > chr22.leg$pos
    nsnps     = sum(gene.snps)

    # status update on parsing
    # only parse genes that span at least 1 SNP
    if (nsnps > 0) {
        cat("Parsing ", sum(gene.snps), " SNPs for gene ", genename, "...\n")
    } else {
        cat("No SNPs found for gene ", genename, ", skipping to next gene...\n")
        next
    }

    # subset genotype files
    aa.geno.sub  = aa.geno[,  gene.snps, with = FALSE]
    ceu.geno.sub = ceu.geno[, gene.snps, with = FALSE]
    yri.geno.sub = yri.geno[, gene.snps, with = FALSE]

    # tack sample information onto the left of the data table 
    # this yields a PLINK RAW format
    aa.geno.sub  = data.table(cbind(aa.samples,  aa.geno.sub))
    ceu.geno.sub = data.table(cbind(ceu.samples, ceu.geno.sub))
    yri.geno.sub = data.table(cbind(yri.samples, yri.geno.sub))

    # add column names for the sample information and the SNPs
    colnames(aa.geno.sub)  = c(colnames(aa.samples),  chr22.leg$rs[gene.snps])
    colnames(ceu.geno.sub) = c(colnames(ceu.samples), chr22.leg$rs[gene.snps])
    colnames(yri.geno.sub) = c(colnames(yri.samples), chr22.leg$rs[gene.snps])

    # make output filepaths for the genotype files
    aa.filename  = file.path(output.dir, paste0("AA.",  genename, ".chr22.raw")) 
    ceu.filename = file.path(output.dir, paste0("CEU.", genename, ".chr22.raw")) 
    yri.filename = file.path(output.dir, paste0("YRI.", genename, ".chr22.raw")) 

    # save subsetted genotypes to file
    fwrite(x = aa.geno.sub,  file = aa.filename,  quote = FALSE, row.names = FALSE, sep = "\t")
    fwrite(x = ceu.geno.sub, file = ceu.filename, quote = FALSE, row.names = FALSE, sep = "\t")
    fwrite(x = yri.geno.sub, file = yri.filename, quote = FALSE, row.names = FALSE, sep = "\t")
}
