#!/usr/bin/env Rscript --vanilla
# ============================================================================================================================== 
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys
#
# This script performs repeated subsampling of GEUVADIS Europeans.
# The resampled number of subjects is always 89, to match the sample size of GEUVADIS Africans.
# ============================================================================================================================== 

suppressMessages(library(methods))
suppressMessages(library(data.table))
suppressMessages(library(optparse))

# error tracking
options(show.error.locations = TRUE)

# parse command line arguments

# parse command line variables
option_list = list(
    make_option(
        c("-p", "--population-name"),
        type    = "character",
        default = NULL,
        help    = "Assign a name (e.g. CEU89) to subsampled population.",
        metavar = "character"
    ),
    make_option(
        c("-e", "--EUR-RNA-file"),
        type    = "character",
        default = NULL,
        help    = "File with RNA data for GEUVADIS Europeans.",
        metavar = "character"
    ),
    make_option(
        c("-a", "--AFR-RNA-file"),
        type    = "character",
        default = NULL,
        help    = "File with RNA data for GEUVADIS Yoruba.",
        metavar = "character"
    ),
    make_option(
        c("-S", "--sample-file"),
        type    = "character",
        default = NULL,
        help    = "PLINK-formatted sample file (two columns, 1 repeated ID per line) of GEUVADIS sample IDs.",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-directory"),
        type    = "character",
        default = NULL,
        help    = "Directory where subsampling output will be saved.",
        metavar = "character"
    ),
    make_option(
        c("-s", "--random-seed"),
        type    = "integer",
        default = 2018,
        help    = "Random seed for reproducible subsamples [default = %default ]",
        metavar = "integer"
    ),
    make_option(
        c("-z", "--subsample-size"),
        type    = "integer",
        default = 89,
        help    = "Size of subsample [default = %default ]",
        metavar = "integer"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

# parse command-line arguments
pop             = opt$population_name    # put a name and sample size to the pop, e.g. YRI89
eur.rna.path    = opt$EUR_RNA_file       # where are the RNA data for the Europeans? 
afr.rna.path    = opt$AFR_RNA_file       # where are the RNA data for the Africans? 
pop.sample.path = opt$sample_file        # where is the PLINK-formatted sample file 
output.dir      = opt$output_directory   # where will pop subsample output go?
seed            = opt$random_seed        # VERY IMPORTANT FOR REPRODUCIBILITY
sample.size     = opt$subsample_size     # make sure that this doesn't exceed the actual sample size 

# set seed
set.seed(seed)  # ALSO VERY IMPORTANT FOR REPRODUCIBILITY

# load expression data
eur.rna = fread(eur.rna.path)
afr.rna = fread(afr.rna.path)

# merge European and African gene expression measurements
rna = merge(eur.rna, afr.rna, by = c("Gene"), all = TRUE)

# load sample IDs for current population
# these are PLINK-formatted files: two columns for IDs
# we only need one of the columns
# nota bene: header = FALSE is VERY IMPORTANT when using fread()
# otherwise the first sample ID becomes the header!
pop.sample.all = fread(pop.sample.path, header = FALSE)[[1]]

# make sure that we only sample from sample IDs actually present in the RNA file
allsamples = c(colnames(eur.rna), colnames(afr.rna))
pop.sample = pop.sample.all[pop.sample.all %in% allsamples] 
cat("length(pop.sample) = ", length(pop.sample), "\n")

# each subsampled set needs its own file paths
pop.rna.path        = file.path(output.dir, paste0("geuvadis.", pop, sample.size, ".RPKM.invnorm.txt"))
pop.pheno.path      = file.path(output.dir, paste0("geuvadis.", pop, sample.size, ".RPKM.invnorm.pheno"))
pop.subjectids.path = file.path(output.dir, paste0("geuvadis.", pop, sample.size, ".sampleids.txt"))

# this variable will house the "Gene" identifiers and a smattering of sample column names
# can use this to (consistently) pull data from eur373.rna
# will also use this to construct files needed for PLINK

new.colnames = c("Gene", sort(pop.sample[sample.int(length(pop.sample), sample.size)]))

# subset the data
pop.rna = rna[, ..new.colnames]

# write this object to file; this is the TXT variant
fwrite(x = pop.rna, file = pop.rna.path, sep = "\t", quote = FALSE)

# melt and recast the data into a PLINK PHENO format
# PHENO format requires FID, IID columns as 1st,2nd columns of data frame
# here, both columns are SubjectID
pop.rna.melt = melt(pop.rna, id.vars = "Gene", variable.name = "SubjectID", value.name = "Measured_Expr")
pop.pheno    = dcast(pop.rna.melt, SubjectID ~ Gene, value.var = "Measured_Expr")
pop.pheno    = cbind(pop.pheno$SubjectID, pop.pheno)

# write PHENO variant to file
# part of PHENO specification is that file DOES NOT have column names
fwrite(x = pop.pheno, file = pop.pheno.path, sep = "\t", quote = FALSE, col.names = FALSE)

# finally, write subject IDs in PLINK-readable format
# this is merely SubjectIDs repeated twice
# as with PHENO, file DOES NOT have column names
fwrite(x = data.table("A" = pop.pheno$SubjectID, "B" = pop.pheno$SubjectID), file = pop.subjectids.path, sep = "\t", quote = FALSE, col.names = FALSE)

# now we repeat this process, but for all IDs *NOT* in the population
# note that we will explicitly exclude members from the population that were not subsampled
# e.g. for CEU, there are 92 samples, of which 89 are subsampled; the remaining 3 are *NOT* included below in the "notceu" group
# but all TSI + GBR + YRI + FIN *ARE* included in "notceu"
notpop.sample          = setdiff(allsamples, pop.sample.all)
notpop.rna.path        = file.path(output.dir, paste0("geuvadis.not", pop, ".RPKM.invnorm.txt"))
notpop.pheno.path      = file.path(output.dir, paste0("geuvadis.not", pop, ".RPKM.invnorm.pheno"))
notpop.subjectids.path = file.path(output.dir, paste0("geuvadis.not", pop, ".sampleids.txt"))

notpop.rna = rna[, ..notpop.sample]
fwrite(x = data.table(notpop.rna), file = notpop.rna.path, sep = "\t", quote = FALSE)

notpop.rna.melt = melt(notpop.rna, id.vars = "Gene", variable.name = "SubjectID", value.name = "Measured_Expr")
notpop.pheno    = dcast(notpop.rna.melt, SubjectID ~ Gene, value.var = "Measured_Expr")
notpop.pheno    = data.table(cbind(notpop.pheno$SubjectID, notpop.pheno))
fwrite(x = notpop.pheno, file = notpop.pheno.path, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(x = data.table("A" = notpop.pheno$SubjectID, "B" = notpop.pheno$SubjectID), file = notpop.subjectids.path, sep = "\t", quote = FALSE, col.names = FALSE)
