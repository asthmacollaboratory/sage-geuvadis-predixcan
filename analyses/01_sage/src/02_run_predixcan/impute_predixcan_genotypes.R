#!/usr/bin/env Rscript --vanilla
#
# this script imputes missing genotype dosaages from a Predixcan genotype file
# it expects two command line arguments:
# -- infile is the original dosage file
# -- outfile is the output path for the new dosage path
#
# coded by Kevin L. Keys (2018)

# parse command line arguments
args    = commandArgs(trailingOnly = TRUE)
infile  = args[1]
outfile = args[2]

# read gzipped file
x = read.table(gzfile(infile), header = FALSE)

# pull out just the dosages (cols 7+)
dosages = x[,-c(1:6)]

# loop over columns of dosages data.frame
# impute missing values of each column to 2*MAF
# put imputed values back into dosages data.frame
for (i in 1:dim(dosages)[2]){
    z = dosages[,i]
    z[is.na(z)] = 2*mean(z, na.rm = T)
    dosages[,i] = z
}

# stick imputed dosages back into data frame
x[,-c(1:6)] = dosages

# write imputed dosages to new gzipped file
write.table(x, file = gzfile(outfile), row.names = FALSE, col.names = FALSE, quote = FALSE)
