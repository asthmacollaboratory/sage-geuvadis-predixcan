#!/usr/bin/env Rscript --vanilla
#
# this little script takes a file of raw RNA-Seq mapped reads and normalizes them to sequencing depth
# the assumed input file structure is a 4-column file with the raw mapped reads in column 4 starting on line 5
# the reads-per-million (RPM) calculation is taken from here
#     http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
# read the file and sum the total mapped reads
# divide that number by 1,000,000 to obtain a reads-per-million scale
# divide the raw reads by that scale to obtain RPMs
# note that GTEx discards samples with < 10M mapped reads, while we do not;
# we need every sample possible, but should at least know about mapping issues

library(data.table)
args         = commandArgs(trailingOnly = TRUE)
infile       = args[1]
x            = fread(infile)
#mapped.reads = x$V4[-c(1:4)]
mapped.reads = x$V4[-c(1,2,4)]
gene.reads   = x$V4[-c(1:4)]
seq.depth    = sum(mapped.reads)
if (seq.depth < 1e7) { warning(paste(infile, "has < 10M mapped reads!\n")) }
norm.scale   = seq.depth / 1e6
rpm          = gene.reads / norm.scale
cat(rpm, "\n")
