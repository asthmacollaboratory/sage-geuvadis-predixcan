##!/usr/bin/env Rscript --vanilla
# ======================================================================================================================
# coded by Kevin L. Keys (2018)
#
# This script downloads a list of genes with start/end positions for a single chromosome.
# It then pads the start/end positions for use in subsetting SNPs in cis regions around the genes. 
# It queries the Ensembl server using the biomaRt package.
# The default is to download genes from chr22 and pad gene boundaries with 500,000 bp. 
# ======================================================================================================================

# ======================================================================================================================
# load libraries 
# ======================================================================================================================
suppressMessages(library(biomaRt))
suppressMessages(library(knitr))
suppressMessages(library(annotate))
suppressMessages(library(data.table))
suppressMessages(library(optparse))

# ======================================================================================================================
# script options
# ======================================================================================================================

# setup, cache = F, echo = FALSE
knitr::opts_chunk$set(error = TRUE)

# ignore warnings, only print errors
options(warn = -1)

# ======================================================================================================================
# command-line arguments 
# ======================================================================================================================

option_list = list(
	make_option(
        c("-o", "--out"),
        type    = "character",
        default = NULL,
        help    = "output file name",
        metavar = "character"
    ),
	make_option(
        c("-a", "--cis-add"),
        type    = "integer",
        default = "500000",
        help    = "cis-region around gene to add to output file, in base pairs [default= %default]",
        metavar = "integer"
    ),
	make_option(
        c("-c", "--chr"),
        type    = "integer",
        default = "22",
        help    = "chromosome from which to grab genes [default= %default]",
        metavar = "integer"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

# ======================================================================================================================
# extract Ensembl IDs, start/end positions for chr22 genes 
# ======================================================================================================================

# get the ensembl mart for humans
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# get list of genes on chr22
gene.list = getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position"),
    filters    = c('chromosome_name'),
    values     = list(opt$chr),
    mart       = ensembl
)

# cast as data.table for easy removal of missing values
gene.list = as.data.table(gene.list)

# set all missing HGNC symbols to NA, then purge rows with NA
gene.list$hgnc_symbol[gene.list$hgnc_symbol == ""] = NA
gene.list = na.omit(gene.list)

# compute lengths of genes
gene.list$gene_length = gene.list$end_position - gene.list$start_position

# add extra region up/downstream of genes for training models
cis.add = opt$cis_add 
gene.list$cis_start = pmax(gene.list$start_position - cis.add, 0)
gene.list$cis_end   = pmin(gene.list$end_position + cis.add, 50818468) ## chr22 is ~50.8Mb long: https://en.wikipedia.org/wiki/Chromosome_22
#gene.list$cis_end   = pmin(gene.list$end_position + cis.add, 248956422) ## longest chr, chr1, is ~249Mb long: https://en.wikipedia.org/wiki/Chromosome_1

# sort gene.list to put longest genes first
setorderv(gene.list, "gene_length", -1L)

# write table to file
fwrite(x = gene.list, file = opt$out, sep = "\t", quote = FALSE)
