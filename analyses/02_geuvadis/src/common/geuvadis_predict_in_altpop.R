##!/usr/bin/env Rscript --vanilla

suppressMessages(library(data.table))
suppressMessages(library(dplyr))

# parse command line arguments
option_list = list(
    make_option(
        c("-b", "--beta-file"),
        type    = "character",
        default = NULL,
        help    = "File with prediction weights (betas) to use in predicting into the test pop.",
        metavar = "character"
    ),
    make_option(
        c("-d", "--genotype-dosage-file"),
        type    = "character",
        default = NULL,
        help    = "File with genotype dosages for the test pop.",
        metavar = "character"
    ),
    make_option(
        c("-g", "--gene-name"),
        type    = "character",
        default = NULL,
        help    = "Name of gene being analyzed.",
        metavar = "character"
    ),
    make_option(
        c("-p", "--prediction-output"),
        type    = "character",
        default = NULL,
        help    = "File where predictions will be saved.",
        metavar = "integer"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

weightsfile      = opt$beta_file 
geno.altpop.file = opt$genotype_dosage_file 
altpop.pred.file = opt$prediction_output 
gene             = opt$gene_name 

# read prediction weights from file
weights = fread(weightsfile)

# purge missing values
weights = na.omit(weights)

# average across all weights for each SNP
betas = data.table(weights %>% group_by(SNP) %>% summarize(Beta = mean(Beta)))

# load PLINK RAW file from alt pop
geno.altpop = fread(geno.altpop.file)

# PLINK RAW normally adds minor/counting allele to marker names
# we need to purge it
colnames(geno.altpop) = gsub("_.*", "", colnames(geno.altpop))

# subset the genotype and weights data.tables to include just the relevant markers
geno.altpop.sub = subset(geno.altpop, select = c("IID", intersect(colnames(geno.altpop), betas$SNP)))
betas.sub       = betas[betas$SNP %in% colnames(geno.altpop.sub),]

# set column order of geno.altpop to match the *row* order of the SNP column in betas.sub 
setorder(betas.sub, SNP)
setcolorder(geno.altpop.sub, c("IID", betas.sub$SNP))

# predictions are now a simple matrix-vector operation
geno.altpop.sub$Prediction_into_Altpop = as.matrix(geno.altpop.sub[,-1]) %*% as.matrix(betas.sub$Beta)

# subset predictions and save to file
altpop.pred      = geno.altpop.sub[,.(IID, Prediction_into_Altpop)]
altpop.pred$Gene = gene
altpop.pred      = altpop.pred[,.(Gene, IID, Prediction_into_Altpop)]
fwrite(x = altpop.pred, file = altpop.pred.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t") 
