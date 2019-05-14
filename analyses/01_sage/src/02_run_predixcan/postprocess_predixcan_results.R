#!/usr/bin/env Rscript --vanilla
#
# this script postprocesses PrediXcan prediction output
# the output is split by chromosome
# unfortunately PrediXcan provides predictions for genes not even on the chromosome
# the script will parse genes that have nonzero predictions and accumulate them into one data frame
# the output is saved based on the inpute file path
#
# this script expects exactly one command line argument:
# -- infile is the prefix to the PrediXcan output files, e.g. "${HOME}/outdir/predixcan_output_" 
#
# coded by Kevin L. Keys (2018)

library(data.table)

# parse command line arguments
args    = commandArgs(trailingOnly = TRUE)
infile  = args[1]

# outfile path is based on input path
outfile = paste(infile, "_ALLCHR_predicted_expression.txt", sep = "")
chrfile = paste(infile, "_chr1_predicted_expression.txt", sep = "")

# alternative file path for melted DF
outfile.melt = paste(infile, "_ALLCHR_predicted_expression_melted.txt", sep = "")

# read first file 
x.chr = fread(chrfile) 

# pull out just the predicted genes 
x.chr.subset = x.chr[,apply(x.chr, 2, function(z) sum(z != 0) != 0), with=FALSE]

# copy chr1 data
# set columns as keys for merging
x = x.chr.subset
setkey(x, FID)
setkey(x, IID)


# trim ENSG IDs to remove transcript number
colnames(x) = strtrim(colnames(x), 15)

# now repeat for remaining chromosomes
for (i in 2:22){
    chrfile = paste(infile, "_chr", i, "_predicted_expression.txt", sep = "")
    x.chr = fread(chrfile) 
    x.chr.subset = x.chr[,apply(x.chr, 2, function(z) sum(z != 0) != 0), with=FALSE]
    setkey(x.chr.subset, FID)
    setkey(x.chr.subset, IID)
    colnames(x.chr.subset) = strtrim(colnames(x.chr.subset), 15)

    # merge current chromosome into accumulated data.table
    x = merge(x, x.chr.subset, by = c("FID", "IID"))
}

# toss FID column 
x = x[, -"FID"]

# now melt the data frame
# rename the columns to something useful
x.melt = melt(x, id.vars = "IID")
colnames(x.melt) = c("IID", "Gene", "Measured_Expr")

# write imputed expressions to file 
fwrite(x=x, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE)
fwrite(x=x.melt, file = outfile.melt, row.names = FALSE, col.names = TRUE, quote = FALSE)
