#!/usr/bin/env Rscript --vanilla
#
# this script subsets RNA-Seq data from gEUVADIS for use with training in glmnet
# the data are pre-downloaded from https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/
# the output of this script consists of six files, three for Europeans (EUR) and two for Yorubans (YRI)
# for each ethnicity, the script produces:
#     (1) a matrix of inverse-normal transformed rank-normalized RPKMs
#     (2) a transposed copy of the matrix from (1)
#     (3) a list of sample IDs for matrices (1) and (2)
#
# This script does not transform the genes to normality. Per the GEUVADIS README
#
#     November 5, 2013 update: The file GD462.GeneQuantRPKM.50FN.samplename.resk10.norm.txt.gz that had the normalization as above
#     PLUS an additional transformation of each gene's values to standard normal has been replaced by GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
#
# https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-3/GeuvadisRNASeqAnalysisFiles_README.txt
#
# coded by Kevin L. Keys (2018)

suppressMessages(library(data.table))
suppressMessages(library(methods))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))

# parse command line variables
option_list = list(
    make_option(
        c("-r", "--rnaseq-data-file"),
        type    = "character",
        default = NULL,
        help    = "File with normalized GEUVADIS gene expression data.",
        metavar = "character"
    ),
    make_option(
        c("-es", "--EUR-sampleIDs-file"),
        type    = "character",
        default = NULL,
        help    = "PLINK-format sample IDs file (two columns, one ID per line, repeated in each column).",
        metavar = "character"
    ),
    make_option(
        c("-ys", "--YRI-sampleIDs-file"),
        type    = "character",
        default = NULL,
        help    = "PLINK-format sample IDs file (two headerless columns, one ID per line, repeated in each column).",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-directory"),
        type    = "character",
        default = NULL,
        help    = "Directory for saving output files",
        metavar = "character"
    ),
    make_option(
        c("-eo", "--EUR-rnaseq-out-file"),
        type    = "character",
        default = "geuvadis.eur373.RPKM.invnorm.txt",
        help    = "Output file for parsed EUR RNA-Seq data [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-yo", "--YRI-rnaseq-out-file"),
        type    = "character",
        default = "geuvadis.yri89.RPKM.invnorm.txt",
        help    = "Output file for parsed YRI RNA-Seq data [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-teo", "--transposed-EUR-rnaseq-out-file"),
        type    = "character",
        default = "geuvadis.eur373.RPKM.invnorm.transposed.txt",
        help    = "Transposed output file for parsed EUR RNA-Seq data [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-tyo", "--transposed-YRI-rnaseq-out-file"),
        type    = "character",
        default = "geuvadis.yri89.RPKM.invnorm.transposed.txt",
        help    = "Transposed output file for parsed YRI RNA-Seq data [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-ei", "--EUR-IDs-out-file"),
        type    = "character",
        default = "geuvadis.eur373.ids.txt",
        help    = "Output file for parsed EUR sample IDs [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-yi", "--YRI-IDs-out-file"),
        type    = "character",
        default = "geuvadis.yri89.ids.txt",
        help    = "Output file for parsed YRI sample IDs [default = %default].",
        metavar = "character"
    )
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

# input file paths
rnaseq.file       = opt$rnaseq_data_file
sample.ids.file   = opt$sample_ids_file
eur.sampleid.file = opt$EUR_sampleIDs_file
yri.sampleid.file = opt$YRI_sampleIDs_file
eur.out.file      = file.path(opt$output_directory, opt$EUR_rnaseq_out_file)
yri.out.file      = file.path(opt$output_directory, opt$YRI_rnaseq_out_file)
teur.out.file     = file.path(opt$output_directory, opt$transposed_EUR_rnaseq_out_file)
tyri.out.file     = file.path(opt$output_directory, opt$transposed_YRI_rnaseq_out_file)
eur.ids.file      = file.path(opt$output_directory, opt$EUR_IDs_out_file)
yri.ids.file      = file.path(opt$output_directory, opt$YRI_IDs_out_file)

# load from file
rnaseq.data     = fread(rnaseq.file, header = TRUE)              # load normalized RPKMs for gEUVADIS
eur.sample.ids  = fread(eur.sampleid.file, header = FALSE)[[1]]  # IDs for EUR, we want these as a *vector*
yri.sample.ids  = fread(yri.sampleid.file, header = FALSE)[[1]]  # same deal for YRI

# discard duplicate gene ID and chr, position info
# rename gene ID to just "Gene"
# also trim any transcript number
rnaseq.data = rnaseq.data[,-c(2:4)]
colnames(rnaseq.data)[1] = "Gene"
rnaseq.data$Gene = strtrim(rnaseq.data$Gene, 15)

# melt data frame
rnaseq.melt = melt(rnaseq.data, value.name = "Expression", variable.name = "IID", id.vars = "Gene")

# parse EUR, YRI from (melted) expression data
# entails filtering IID column against relevant sample IDs
# recast result into wide format
eur = rnaseq.melt %>%
    dplyr::filter(IID %in% eur.sample.ids) %>%
    dcast(Gene ~ IID, value.var = "Expression") %>%
    as.data.table

yri = rnaseq.melt %>%
    dplyr::filter(IID %in% yri.sample.ids) %>%
    dcast(Gene ~ IID, value.var = "Expression") %>%
    as.data.table

# order eur and yri by gene ID
setorder(x=yri, Gene, na.last = TRUE)
setorder(x=eur, Gene, na.last = TRUE)

# write subsetted data.tables to file
fwrite(x = eur, file = eur.out.file, col.names = TRUE, quote = FALSE, sep = "\t")
fwrite(x = yri, file = yri.out.file, col.names = TRUE, quote = FALSE, sep = "\t")

# also write transposed data files
# involves transposing data tables via melting and recasting
# do first for eur and then for yri
eur.melt = melt(eur, id.vars = "Gene", variable.name = "IID", value.name = "Expression")
teur     = dcast(data = eur.melt, formula = IID ~ Gene, value.var = "Expression")
setkey(teur, IID)
setorder(teur, IID)

yri.melt = melt(yri, id.vars = "Gene", variable.name = "IID", value.name = "Expression")
tyri     = dcast(data = yri.melt, formula = IID ~ Gene, value.var = "Expression")
setkey(tyri, IID)
setorder(tyri, IID)

fwrite(x = teur, file = teur.out.file, quote = FALSE, sep = "\t")
fwrite(x = tyri, file = tyri.out.file, quote = FALSE, sep = "\t")

# also save separate lists of EUR, YRI sample IDs
# this will come in handy when subsetting genotypes
fwrite(x = data.table(eur.sample.ids), file = eur.ids.file, col.names = FALSE, quote = FALSE, sep = "\t")
fwrite(x = data.table(yri.sample.ids), file = yri.ids.file, col.names = FALSE, quote = FALSE, sep = "\t")
