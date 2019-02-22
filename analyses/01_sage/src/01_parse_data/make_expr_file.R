#!/usr/bin/env Rscript --vanilla
#
# this script operates on the RNA-Seq read counts from 39 subjects from SAGE
# it normalizes the *sequence-depth-normalized read counts* to RPKM and outputs the normalized expression levels
# the script is meant to be used from the command line ONLY and requires four arguments in the following order:
#   1) a file path to the input matrix of RPMs 
#   2) a file path to save the expression file
#   3) a file path to save the PLINK-formatted ID file
#   4) a file path to the input matrix of gene information, which includes gene length 
#   5) a file path to save the RPKMs computed from the RPMs in (1)
#
# an example call:
#
# > Rscript make_matrixeqtl_expr.file.R $READS $OUTPUT $IDS $GENES
#
# coded by Kevin L. Keys (2017)

# load libraries
library(methods)
library(data.table)
library(preprocessCore)
library(optparse)

# script variables used for data cleaning
# these numbers are drawn from the GTEx pipeline v6p
# -- > 0.1 RPKM in at least 10 individuals
# -- > 5 reads in at least 10 individuals
rpkm.threshold  = 0.1
count.threshold = 5
#min.samples     = 10 # GTEx v6: normalize based on fixed number of samples
min.samples     = 0.2*39 # GTEx v7: normalized based on percentage of samples

### new 6 FEB 2018:
# after speaking with Walter Eckalbar, we will relax these thresholding criteria
# our sample size is small enough to make these thresholds problematic (they discard too many data points)
rpkm.threshold = 1e-6
count.threshold = 5
min.samples = 0.2*39

# parse command line arguments
option_list = list(
    make_option(
        c("-a", "--RPM-file"),
        type    = "character",
        default = NULL,
        help    = "File of reads per million",
        metavar = "character"
    ),
    make_option(
        c("-b", "--output-file"),
        type    = "character",
        default = NULL,
        help    = "Output file for normalized expression data",
        metavar = "character"
    ),
    make_option(
        c("-c", "--ID-outfile"),
        type    = "character",
        default = NULL,
        help    = "Output file for IDs",
        metavar = "character"
    ),
    make_option(
        c("-d", "--genefile"),
        type    = "character",
        default = NULL,
        help    = "List of genes downloaded from Ensembl",
        metavar = "character"
    ),
    make_option(
        c("-e", "--RPKM-file"),
        type    = "character",
        default = NULL,
        help    = "Output file of RPKM expression values",
        metavar = "character"
    ),
    make_option(
        c("-f", "--reads-file"),
        type    = "character",
        default = NULL,
        help    = "File containing RNA-Seq reads",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

# file paths
my.file      = opt$RPM_file    # sage_rnaseq_rpms.txt
my.outfile   = opt$output_file 
my.idfile    = opt$ID_outfile 
my.genefile  = opt$genefile    # human_ens_GRCh37_annot.extended.txt
my.rpkmfile  = opt$RPKM_file 
my.readsfile = opt$reads_file  # sage_rnaseq_reads.txt 

# load data
x     = fread(my.file)
y     = fread(my.genefile)
reeds = fread(my.readsfile)

# transpose data
xt = data.frame(t(x[,-1]))
colnames(xt) = gsub("out", "", x$SubjectID)
#x$Ensembl_ID = row.names(xt)
treeds = data.frame(t(reeds[,-1]))
colnames(treeds) = colnames(xt)

# subset y based on rows of xt
y.subset = subset(y, y$Ensembl_ID %in% row.names(xt))
setorder(y.subset, Ensembl_ID)

# normalize xt based on gene lengths
xt = xt[order(row.names(xt)),]
if(sum(row.names(xt) != y.subset$Ensembl_ID) > 0){ warning("unmatched total exon lengths and gene IDs when normalizing expression data!") }
xt.norm = xt / y.subset$Total_Exon_Length
treeds = treeds[order(row.names(treeds)),]

# make an RPKMs data frame to save
rpkms = xt.norm
rpkms$Gene = y.subset$Ensembl_ID
rpkms = rpkms[,c(ncol(rpkms),1:(ncol(rpkms)-1))]

# write RPKMs to file
write.table(
    file  = my.rpkmfile,
    x     = rpkms, 
    sep   = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE 
)

# subroutine to compute rank
# this subroutine does not count NA in the ranking
# it ensures that NA remains NA after ranking
# see https://stackoverflow.com/questions/19409550/rank-and-length-with-missing-values-in-r
# update: the function
# > prank = function(z) ifelse(is.na(z), NA, rank(z) / sum(!is.na(z)) )
# is provided in R as
# ecdf(x)(x)
# see Hadley Wickham's comment at aforementioned SO webpage
prank = function(z) ecdf(z)(z)

# select expression levels based on given thresholds 
mask     = (apply(xt.norm > rpkm.threshold, 1, sum) >= min.samples) & (apply(treeds > count.threshold, 1, sum) >= min.samples)
temp.x   = xt.norm[mask,]
temp.x2  = normalize.quantiles(as.matrix(temp.x)) # output of quantile normalization: rows = genes, columns = samples
temp.x3  = apply(temp.x2, 1, prank) # output of rank computation: rows = SAMPLES, cols = GENES
xt.norm2 = data.frame(qnorm(t(temp.x3))) # qnorm operates on columns. equivalently, xt.norm2 = data.frame(apply(temp.x3, 2, qnorm))
colnames(xt.norm2) = colnames(xt.norm)

# add column of "Gene" names to xt.norm
# this ensures readable headers later
xt.norm3 = cbind("Gene" = row.names(temp.x), xt.norm2)

# write transposed data to file
write.table(
    file  = my.outfile,
    x     = xt.norm3,
    sep   = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

# also want to save subject IDs for subsetting with PLINK
write.table(
    file  = my.idfile,
    x     = cbind(row.names(x), row.names(x)),
    quote = F,
    row.names = FALSE,
    col.names = FALSE
)
