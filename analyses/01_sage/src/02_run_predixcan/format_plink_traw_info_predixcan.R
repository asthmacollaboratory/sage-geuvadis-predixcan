#!/usr/bin/env Rscript --vanilla
# This script formats input file for PrediXcan
# coded by Kevin L. Keys (2018)
# based on script by Donglei Hu

library(methods)

# main subroutine for this script
# it loads a PLINK TRAW and a FREQ file
# and then creates a PrediXcan-compatible input file 
format.plink.traw.into.predixcan = function(inNameGeno, inNameFreq, outName){

    # load files
    geno = read.table(inNameGeno,header=T,stringsAsFactors=F)
    freq = read.table(inNameFreq,header=T,stringsAsFactors=F)

    # check that the number of SNPs in each file is the same
    # also ensure that the number of counted minor alleles is the same in both files 
    if ( (sum(geno$SNP==freq$SNP)==nrow(geno)) & (sum(geno$COUNTED==freq$A1)==nrow(freq)) ) {

        # parse the current chromosome
        chr = paste('chr', freq$CHR, sep='')

        # here we insert the MAF before the dosages in the TRAW file
        predixcan = cbind(chr, geno[,c(2,4,6,5)], freq[,5], geno[,7:ncol(geno)])

        # write the output to file
        # pass it through gzip to save space
        write.table(x = predixcan, file = gzfile(outName), quote=FALSE, row.names=FALSE, col.names=FALSE)
    }

    return()
}

# parse command line arguments
# There are 2 input files
# Input file 1: plink traw file
# Input file 2: plink frq file
# the third argument specifies the output path
args = commandArgs(trailingONLY = TRUE)
infile.geno = args[1]
infile.freq = args[2]
outfile     = args[3]
    
# execute reformatting
format.plink.traw.into.predixcan()

#inNameGeno='Geno_bdr_1441_refiltered_rsid_tw_wb_chr19.traw'
#inNameFreq='Freq_bdr_1441_refiltered_rsid_tw_wb_chr19.frq'
#outName='chr19_bdr_1441_refiltered_tw_wb_predixcan.gz'
