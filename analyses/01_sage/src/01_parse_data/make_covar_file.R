#!/usr/bin/env Rscript --vanilla
#
# this script constructs the covariate file for MatrixEQTL
#
# example:
#
# > Rscript make_matrixeqtl_covar_file.R $IDFILE $PHENOFILE $COVARFILE
#
# coded by Kevin L. Keys (2017)

# script variables
covars.to.include = c("SubjectID", "age", "Sex", "bmi")
new.covar.names   = c("SubjectID", "Age", "Sex", "BMI")

# load libraries
library(methods)
library(data.table)
library(optparse)

# parse command line arguments
option_list = list(
    make_option(
        c("-a", "--phenotype-file"),
        type    = "character",
        default = NULL,
        help    = "File of sample phenotypes.",
        metavar = "character"
    ),
    make_option(
        c("-b", "--output-file"),
        type    = "character",
        default = NULL,
        help    = "Output file for covariates.",
        metavar = "character"
    ),
    make_option(
        c("-d", "--PC-file"),
        type    = "character",
        default = NULL,
        help    = "File with principal components from PLINK.",
        metavar = "character"
    ),
    make_option(
        c("-e", "--PEER-file"),
        type    = "character",
        default = NULL,
        help    = "Output file from R PEER computations.",
        metavar = "character"
    ),
    make_option(
        c("-f", "--transposed-output-file"),
        type    = "character",
        default = NULL,
        help    = "Output file for *transposed* covariates.",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

# parse command line arguments
phenotypes     = opt$phenotype_file
pcfile         = opt$PC_file
outfile        = opt$output_file
peerfile       = opt$PEER_file 
toutfile       = opt$transposed_output_file

# load phenotype file
x = fread(phenotypes, header = TRUE)

# load PEER covariates; the GTEx pipeline uses a variable number of these
# for N = 39 in our case, we need fifteen PEER factors
# see https://gtexportal.org/home/documentationPage for more info
# this file is output from the PEER package in R and has a header and row names
# the row names are terrible so rename them 
peer = read.table(peerfile, header=TRUE, row.names=1)
for (i in 1:dim(peer)[1]){ 
    rownames(peer)[i] = paste("PEER", i, sep = "")
}


# load PC file
# add names to the columns
# this assumes that we have only 3 PCs!
pc.file = fread(pcfile, header = FALSE)
colnames(pc.file) = c("SubjectID", "IID", "PC1", "PC2", "PC3")
pc.file = subset(pc.file, select = c("SubjectID", "PC1", "PC2", "PC3"))

# subset the covariate file
covars = subset(x, idx, select = covars.to.include)
colnames(covars) = new.covar.names

# merge the covariates with the PCs 
covars = merge(covars, pc.file, by = "SubjectID") 

# need to enforce numeric dummy variables
# this only applies to sex
# with this encoding, 0 --> FEMALE and 1 --> MALE
covars$Sex = as.numeric(as.factor(covars$Sex)) - 1

# transpose the file without the SubjectIDs
tcovars = t(as.data.frame(covars[,-1]))
colnames(tcovars) = covars$SubjectID

# lastly, slap the PEER factors onto the bottom of the covariate data frame
tcovars = rbind(tcovars, peer)

# write covariate file to file 
write.table(
    x = tcovars,
    file = outfile,
    quote = FALSE,
    row.names = TRUE,
    sep = "\t"
)

tpeer = as.data.frame(t(peer))
tpeer$SubjectID = row.names(tpeer)
covars = merge(covars, tpeer, by = "SubjectID")
colnames(covars)[1] = "IID"


# write "transposed" covariate file to file 
write.table(
    x = covars,
    file = toutfile,
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
)
