##!/usr/bin/env Rscript --vanilla

suppressMessages(library(data.table))
suppressMessages(library(methods))
suppressMessages(library(optparse))

# parse command line variables
option_list = list(
    make_option(
        c("-f", "--bim-file"),
        type    = "character",
        default = NULL,
        help    = "PLINK BIM file for GEUVADIS genotype data",
        metavar = "character"
    ),
    make_option(
        c("-i", "--SNP-ID-file"),
        type    = "character",
        default = NULL,
        help    = "List of dbSNP IDs to use for remapping",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-file"),
        type    = "character",
        default = NULL,
        help    = "Path to file for saving results.",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

# parse command line arguments
in.file  = opt$bim_file
id.file  = opt$SNP_ID_file
out.file = opt$output_file

# load PLINK BIM
x = fread(in.file, header = FALSE)

# load ID conversion file
y = fread(id.file, header = FALSE)

# remap gEUVADIS SNP IDs to rsIDs
rsids = y$V1[y$V2 %in% x$V2]
x$V2  = rsids

# write new BIM to file
fwrite(x = x, file = out.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
