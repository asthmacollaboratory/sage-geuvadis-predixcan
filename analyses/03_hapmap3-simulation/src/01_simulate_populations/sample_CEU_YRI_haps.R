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
    ),
    make_option(
        c("-z1", "--num-CEU"),
        type    = "integer",
        default = 10000, 
        help    = "Number of samples to generate from CEU haplotypes [default = %default]",
        metavar = "integer"
    ),
    make_option(
        c("-z2", "--num-YRI"),
        type    = "integer",
        default = 10000, 
        help    = "Number of samples to generate from YRI haplotypes [default = %default]",
        metavar = "integer"
    ),
    make_option(
        c("-z3", "--num-AA"),
        type    = "integer",
        default = 10000, 
        help    = "Number of AA samples to generate from CEU/YRI haplotypes [default = %default]",
        metavar = "integer"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

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
pad               = as.integer(opt$pad_around_genes)
ngenes            = as.integer(opt$num_genes)
ceu.num           = as.integer(opt$num_CEU)
yri.num           = as.integer(opt$num_YRI)
aa.num            = as.integer(opt$num_AA)

# set random seed
set.seed(seed)

# set output file paths
# these geno files will contain all of the allele dosages for all chr22 genes
aa.geno.file  = file.path(output.dir, paste0("AA.chr22.CEU_", ceu.prop, ".YRI_", yri.prop, ".geno"))
ceu.geno.file = file.path(output.dir, "CEU.chr22.geno")
yri.geno.file = file.path(output.dir, "YRI.chr22.geno")

# load haps
ceu.haps = fread(ceu.hap.file.path)
yri.haps = fread(yri.hap.file.path)

## sample proportions
#ceu.num = dim(ceu.haps)[2]
#yri.num = dim(yri.haps)[2]
#ceu.haps.sub = ceu.haps[,sample(1:ceu.num,floor(ceu.prop*ceu.num)), with = FALSE]
#yri.haps.sub = yri.haps[,sample(1:yri.num,ceiling(yri.prop*yri.num)), with = FALSE]

# sample proportions
sample.ceu.replacement = FALSE
sample.yri.replacement = FALSE
sample.aa.replacement  = FALSE
ceu.num.haps = ncol(ceu.haps)
yri.num.haps = ncol(yri.haps)

# this allows for sampling a # of CEU haps that differs from what is in the hap file
# everything is copesetic if the user asks for the # of CEU haps in the hap file
# if user asks for more, then we will sample with replacement
# in contrary case, will simply sample fewer haps
if (ceu.num != 2*ceu.num.haps) {
    if (ceu.num > 2*ceu.num.haps) {
        sample.ceu.replacement = TRUE
    }
    #ceu.num.haps = 2*as.numeric(opt$num_CEU)
    ceu.num.haps = min(ceu.num.haps, 2*as.numeric(opt$num_CEU))
}

# do same sampling adjustment for YRI
if (yri.num != 2*yri.num.haps) {
    if (yri.num > 2*yri.num.haps) {
        sample.yri.replacement = TRUE
    }
    #yri.num.haps = 2*(opt$num_YRI)
    yri.num.haps = min(yri.num.haps, 2*(opt$num_YRI))
}

# similar sampling adjustment for AA
# this one depends on CEU+YRI haps
aa.num.haps = floor(ceu.prop*ceu.num.haps) + ceiling(yri.prop*yri.num.haps) 
if (aa.num != 2*aa.num.haps) {
    if (aa.num > 2*aa.num.haps) {
        sample.aa.replacement = TRUE
    }
    aa.num.haps = min(aa.num.haps, 2*(opt$num_AA))
}

cat("will generate N = ", ceu.num, " CEU samples from ", ceu.num.haps, "haplotypes\n")
cat("will generate N = ", yri.num, " YRI samples from ", yri.num.haps, "haplotypes\n")
cat("will generate N = ", aa.num,  " AA samples from ",  aa.num.haps,  "haplotypes\n")

# sample the haplotypes
#ceu.haps.sub = ceu.haps[,sample(1:ceu.num.haps, floor(ceu.prop*ceu.num.haps), replace = sample.ceu.replacement), with = FALSE]
#yri.haps.sub = yri.haps[,sample(1:yri.num.haps, ceiling(yri.prop*yri.num.haps), replace = sample.yri.replacement), with = FALSE]
ceu.haps.sample.for.aa = sample.int(ceu.num.haps, size = floor(ceu.prop*ceu.num.haps), replace = sample.ceu.replacement)
yri.haps.sample.for.aa = sample.int(yri.num.haps, size = ceiling(yri.prop*yri.num.haps), replace = sample.yri.replacement)
ceu.haps.sub.for.aa = ceu.haps[,ceu.haps.sample.for.aa, with = FALSE]
yri.haps.sub.for.aa = yri.haps[,yri.haps.sample.for.aa, with = FALSE]
ceu.haps.sample.file.for.aa = file.path(output.dir, paste0("AA.fromCEU.chr22.CEU_", ceu.prop, ".YRI_", yri.prop, ".haps.txt"))
yri.haps.sample.file.for.aa = file.path(output.dir, paste0("AA.fromYRI.chr22.CEU_", ceu.prop, ".YRI_", yri.prop, ".haps.txt"))
fwrite(x = data.table(ceu.haps.sample.for.aa), file = ceu.haps.sample.file.for.aa, quote = FALSE, row.names = FALSE, col.names = FALSE)
fwrite(x = data.table(yri.haps.sample.for.aa), file = yri.haps.sample.file.for.aa, quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("sampled ", ncol(ceu.haps.sub.for.aa), " haplotypes for AA from CEU per proportion ", ceu.prop, "\n")
cat("sampled ", ncol(yri.haps.sub.for.aa), " haplotypes for AA from YRI per proportion ", yri.prop, "\n")

# generate a sample of the AA haps - this will shuffle the CEU/YRI haplotypes that make AA
# without this, then AA would just be $ceu.prop CEU ppl and $yri.prop YRI ppl
aa.haps.sample = sample.int(aa.num.haps, size = aa.num.haps, replace = sample.aa.replacement)

# save the haplotype sample order for bookkeeping
aa.haps.sample.file = file.path(output.dir, paste0("AA.chr22.CEU_", ceu.prop, ".YRI_", yri.prop, ".haps.txt"))
fwrite(x = data.table(aa.haps.sample), file = aa.haps.sample.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

# slap haps together to make AA
if (ceu.prop == 0) {
    aa.haps = data.table(yri.haps.sub.for.aa)[, aa.haps.sample, with = FALSE]
} else if (yri.prop == 0) {
    aa.haps = data.table(ceu.haps.sub.for.aa)[, aa.haps.sample, with = FALSE]
} else {
    ###aa.haps = data.table(cbind(ceu.haps.sub, yri.haps.sub))
    #aa.num.haps = ncol(ceu.haps.sub) + ncol(yri.haps.sub)
    #aa.haps = data.table(cbind(ceu.haps.sub, yri.haps.sub))[,sample.int(aa.num.haps, size = aa.num.haps, replace = FALSE), with = FALSE]

    # slap the haplotypes together and then shuffle
    aa.haps = data.table(cbind(ceu.haps.sub.for.aa, yri.haps.sub.for.aa))[, aa.haps.sample, with = FALSE]
}

# keep track of which haplotype came from where (CEU or YRI)
aa.haps.origin = c(rep("CEU", ncol(ceu.haps.sub.for.aa)), rep("YRI", ncol(yri.haps.sub.for.aa)))[aa.haps.sample]

# save the haplotype origin to file too
aa.haps.origin.file = file.path(output.dir, paste0("AA.chr22.CEU_", ceu.prop, ".YRI_", yri.prop, ".haps.origin.txt"))
fwrite(x = data.table(aa.haps.origin), file = aa.haps.origin.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

aa.num.haps  = ncol(aa.haps)
aa.num = min(ceu.num, yri.num)
cat("will generate N = ", aa.num, "samples of AA drawn from a total of ", aa.num.haps, " haplotypes\n")

# write to file
aa.haps.file = file.path(output.dir, paste0("AA.chr22.CEU_", ceu.prop, ".YRI_", yri.prop, ".haps"))
fwrite(x = aa.haps, file = aa.haps.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

# create genotype files
# do this by collapsing every other column 
ceu.geno = data.table(ceu.haps[,seq(1,ceu.num.haps,2), with = FALSE] + ceu.haps[,seq(2,ceu.num.haps,2), with = FALSE])
yri.geno = data.table(yri.haps[,seq(1,yri.num.haps,2), with = FALSE] + yri.haps[,seq(2,yri.num.haps,2), with = FALSE])
aa.geno  = data.table(aa.haps[,seq(1,aa.num.haps,2), with = FALSE]   + aa.haps[,seq(2,aa.num.haps,2), with = FALSE])
gc()

ceu.haps.sample = sample.int(ceu.num.haps, size = 2*ceu.num, replace = sample.ceu.replacement)
yri.haps.sample = sample.int(yri.num.haps, size = 2*yri.num, replace = sample.yri.replacement)
ceu.haps.sub = ceu.haps[,ceu.haps.sample, with = FALSE]
yri.haps.sub = yri.haps[,yri.haps.sample, with = FALSE]
ceu.haps.sample.file = file.path(output.dir, paste0("CEU.chr22.CEU_", ceu.prop, ".YRI_", yri.prop, ".haps.txt"))
yri.haps.sample.file = file.path(output.dir, paste0("YRI.chr22.CEU_", ceu.prop, ".YRI_", yri.prop, ".haps.txt"))
fwrite(x = data.table(ceu.haps.sample), file = ceu.haps.sample.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
fwrite(x = data.table(yri.haps.sample), file = yri.haps.sample.file, quote = FALSE, row.names = FALSE, col.names = FALSE)

ceu.geno = data.table(ceu.haps.sub[,seq(1,ncol(ceu.haps.sub),2), with = FALSE] + ceu.haps[,seq(2,ncol(ceu.haps.sub),2), with = FALSE])
yri.geno = data.table(yri.haps.sub[,seq(1,ncol(yri.haps.sub),2), with = FALSE] + yri.haps[,seq(2,ncol(yri.haps.sub),2), with = FALSE])
aa.geno  = data.table(aa.haps[,seq(1,ncol(aa.haps),2), with = FALSE]   + aa.haps[,seq(2,ncol(aa.haps),2), with = FALSE])
gc()

# similar to haplotypes, keep track of genotype origins
# this of this file as a crude local ancestry estimator
aa.geno.origin = paste(aa.haps.origin[seq(1,length(aa.haps.origin),2)], aa.haps.origin[seq(2,length(aa.haps.origin),2)], sep = "/") 

# save the haplotype origin to file too
aa.geno.origin.file = file.path(output.dir, paste0("AA.chr22.CEU_", ceu.prop, ".YRI_", yri.prop, ".geno.origin.txt"))
fwrite(x = data.table(aa.geno.origin), file = aa.geno.origin.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

# the $*.geno files contain SNPs in rows and genotypes in columns
# transpose to SNP-major format (SNPs on columns)
ceu.geno = data.table(t(ceu.geno))
yri.geno = data.table(t(yri.geno))
aa.geno  = data.table(t(aa.geno))

cat("made CEU genotype file with ", nrow(ceu.geno), " samples and ", ncol(ceu.geno), " SNPs\n") 
cat("made YRI genotype file with ", nrow(yri.geno), " samples and ", ncol(yri.geno), " SNPs\n") 
cat("made AA  genotype file with ", nrow(aa.geno), " samples and ", ncol(aa.geno), " SNPs\n") 

# write genotype files to file
# only write CEU/YRI genotypes if they do not already exist
# this saves some time by avoiding overwriting the same file with the same information 
fwrite(x = aa.geno,  file = aa.geno.file,  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
if (!file.exists(ceu.geno.file)) {
    fwrite(x = ceu.geno, file = ceu.geno.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}
if (!file.exists(yri.geno.file)) {
    fwrite(x = yri.geno, file = yri.geno.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

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
    aa.filename  = file.path(output.dir, paste0("AA.",  genename, "_CEU", ceu.prop, "_YRI", yri.prop, ".chr22.raw")) 
    ceu.filename = file.path(output.dir, paste0("CEU.", genename, ".chr22.raw")) 
    yri.filename = file.path(output.dir, paste0("YRI.", genename, ".chr22.raw")) 

    # save subsetted genotypes to file
    fwrite(x = aa.geno.sub,  file = aa.filename,  quote = FALSE, row.names = FALSE, sep = "\t")
    if (!file.exists(ceu.filename)) {
        fwrite(x = ceu.geno.sub, file = ceu.filename, quote = FALSE, row.names = FALSE, sep = "\t")
    }
    if (!file.exists(yri.filename)) {
        fwrite(x = yri.geno.sub, file = yri.filename, quote = FALSE, row.names = FALSE, sep = "\t")
    }
}
