#!/usr/bin/env Rscript --vanilla
# ========================================================================================
# coded by Kevin L. Keys
#
# This script subsets populations from the GEUVADIS sample set.
# It outputs PLINK-compatible ID lists for subsetting the GEUVADIS genotype data.
#
# Subsets and corresponding output file names include:
# -- individual populations: "geuvadis.${POP}${POP_SIZE}.sampleids.txt"
# -- complement of individual populations: "geuvadis.not${POP}.sampleids.txt
# -- the EUR278 set excluding Africans and Finnish
#
# For analysis, note the following name changes:
# -- geuvadis.notyri.sampleids.txt ==> geuvadis.eur373.sampleids.txt
#
# ========================================================================================
suppressMessages(library(data.table))

# read data from 1000Genomes website: http://www.internationalgenome.org/data-portal/sample
# data were generated by
# 1. Clicking "Filter by data collection"
# 2. Checking "Geuvadis"
# 3. Clicking "Download the list"
datafiles.dir = "../../datafiles"
path.to.samples = file.path(datafiles.dir, "igsr_samples.tsv")
x = fread(path.to.samples)

# GEUVADIS has 5 constituent populations
pops = c("GBR", "TSI", "CEU", "YRI", "FIN")

# parse the populations themselves
for (pop in pops) {
	x.sub = x[c(pop),c("Sample name", "Sample name")]
	x.num = dim(x.sub)[1]
	filename = filepath(datafiles.dir, paste0("geuvadis.", tolower(pop), x.num, ".sampleids.txt"))
	fwrite(x.sub, file = filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

# parse the complement of the populations
# e.g. for Africans, parse "NOT YRI" a.k.a. the Europeans
for (pop in pops) {
	x.sub = x[!c(pop),c("Sample name", "Sample name")]
	x.num = dim(x.sub)[1]
	filename = file.path(datafiles.dir, paste0("geuvadis.not", tolower(pop), ".sampleids.txt"))
	fwrite(x.sub, file = filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

# parse the complement of YRI + FIN
x.sub = x[!c("YRI", "FIN"),c("Sample name", "Sample name")]
filename = file.path(datafiles.dir, "geuvadis.eur278.sampleids.txt")
fwrite(x.sub, file = filename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
