##!/usr/bin/env Rscript --vanilla

# ==============================================================================================================================
# load libraries
# ==============================================================================================================================
suppressMessages(library(dplyr))
suppressMessages(library(data.table))


# ==============================================================================================================================
# parse command line arguments 
# ==============================================================================================================================

option_list = list(
    make_option(
        c("-b", "--beta-file"),
        type    = "character",
        default = NULL,
        help    = "File with prediction weights (betas).",
        metavar = "character"
    ),
    make_option(
        c("-d", "--discard-ratio"),
        type    = "double",
        default = NULL,
        help    = "Minimum admissible ratio of predictions for analyzing a gene.",
        metavar = "double"
    ),
    make_option(
        c("-np", "--num-predictions-file"),
        type    = "character",
        default = NULL,
        help    = "File to save numbers of predictions computed per gene.", 
        metavar = "character"
    ),
    make_option(
        c("-ns", "--num-samples-file"),
        type    = "character",
        default = NULL,
        help    = "File where predictions from glmnet will be saved.",
        metavar = "character"
    ),
    make_option(
        c("-s", "--num-samples"),
        type    = "character",
        default = NULL,
        help    = "Number of samples in the training set.",
        metavar = "character"
    ),
    make_option(
        c("-p", "--prediction-file"),
        type    = "character",
        default = NULL,
        help    = "File with predictions for the training population.", 
        metavar = "character"
    ),
    make_option(
        c("-r", "--expression-file"),
        type    = "character",
        default = NULL,
        help    = "File with table of gene expression measurements in the training population.",
        metavar = "character"
    ),
    make_option(
        c("-ol", "--out-lm-file"),
        type    = "character",
        default = NULL,
        help    = "Output file for saving linear modeling results across all genes.",
        metavar = "character"
    ),
    make_option(
        c("-og", "--out-genelm-file"),
        type    = "character",
        default = NULL,
        help    = "Output file for saving a table of linear modeling results for each gene.",
        metavar = "character"
    ),
    make_option(
        c("-tp", "--test-pop-prediction-file"),
        type    = "character",
        default = NULL,
        help    = "File with predictions for the testing population.",
        metavar = "character"
    ),
    make_option(
        c("-tr", "--test-pop-expression-file"),
        type    = "character",
        default = NULL,
        help    = "File with gene expression measurements for the testing population.",
        metavar = "character"
    ),
    make_option(
        c("-tol", "--test-pop-out-lm-file"),
        type    = "character",
        default = NULL,
        help    = "Output file for saving linear modeling results across genes in the testing population.",
        metavar = "character"
    ),
    make_option(
        c("-tog", "--test-pop-out-genelm-file"),
        type    = "character",
        default = NULL,
        help    = "Output file for saving a table of linear modeling results for each gene in the testing population.",
        metavar = "character"
    )
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

# parse command line arguments
weightsfile            = opt$beta_file 
#new.weightsfile        = args[2]
discard.ratio          = opt$discard_ratio 
num.pred.file          = opt$num_predictions_file 
num.samples            = opt$num_samples 
prediction.file        = opt$prediction_file
exprfile               = opt$expression_file 
out.lm.file            = opt$out_lm_file 
out.genelm.file        = opt$out_genelm_file 
altpop.pred.file       = opt$test_pop_prediction_file 
altpop.exprfile        = opt$test_pop_expression_file 
altpop.out.lm.file     = opt$test_pop_out_lm_file 
altpop.out.genelm.file = opt$test_pop_out_genelm_file 

# ensure that numeric arguments are actually numbers
discard.ratio = as.numeric(discard.ratio)
num.samples   = as.numeric(num.samples)

# ==============================================================================================================================
# subroutines
# ==============================================================================================================================

compute.gene.lm = function(data, out.file, genes = unique(sort(data$Gene))){
    ngenes = length(genes)
    out.matrix     = matrix(NA, ngenes, 4)
    out.matrix[,1] = genes
    for (i in 1:ngenes) {
        gene = genes[i]
        data.sub = data[data$Gene == gene,]
        if(!all(is.na(data.sub$Predicted_Expr))){
            my.lm  = summary(lm(Predicted_Expr ~ Measured_Expr, data = data.sub))
            lm.p   = my.lm$coefficients[2,4]
            lm.r2  = my.lm$r.squared
            my.rho = cor.test(data.sub$Predicted_Expr, data.sub$Measured_Expr, method = "spearman")$estimate
            out.matrix[i,2:4] = c(lm.p, lm.r2, my.rho)
        }   
    }   
    out.df = data.table(out.matrix)
    colnames(out.df) = c("Gene", "P.value", "R2", "Spearman.rho")
    fwrite(out.df, file = out.file, row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA", sep = "\t")
    return()
}

# ==============================================================================================================================
# script code 
# ==============================================================================================================================

# open weights file
x = fread(weightsfile)

# ensure proper header
my.header = c("Gene","Held_out_Sample","SNP","A1","A2","Beta")
colnames(x) = my.header

# determine how many predicted samples that each gene has
# we will sink this into an output file and reload it as a data.table
cat(paste("creating num.pred.file at ", num.pred.file, "...\n", sep = ""))
num.pred = data.table(count(count(x, Gene, Held_out_Sample), Gene))
colnames(num.pred) = c("Gene", "Num.Pred")
fwrite(x = num.pred, file = num.pred.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t") 
cat("done.\n")

# subset those genes with at least discard.ratio * (# samples) predicted samples
cat("discarding genes with less than ", discard.ratio * num.samples, "predictions...\n")
num.pred.sub = num.pred[num.pred$Num.Pred > discard.ratio*num.samples,]

# we can now look at the predictions for the "well-predicted" genes in num.pred.sub
# must load prediction information first
cat("reading prediction.file at ", prediction.file, "...\n")
geuvadis.predictions = fread(prediction.file)

# now load measured RNA data
cat("reading gene expression file at ", exprfile, "...\n")
geuvadis.rnaseq = fread(exprfile)


# the rnaseq and prediction files have the same column name order
# expected output, using YRI 89 as an example:
# > colnames(geuvadis.rnaseq)
#     [1] "Gene"    "NA18486" "NA18487" "NA18488" "NA18489" "NA18498" "NA18499"
#     [8] ... 
colnames(geuvadis.predictions) = colnames(geuvadis.rnaseq)

# now subset predictions based on which genes are well predicted
geuvadis.predictions.sub = geuvadis.predictions[geuvadis.predictions$Gene %in% num.pred.sub$Gene,]

# subset the file with just the genes from geuvadis.predictions.sub
cat("subsetting expression file to well-predicted genes...\n") 
geuvadis.rnaseq.sub = geuvadis.rnaseq[geuvadis.rnaseq$Gene %in% geuvadis.predictions.sub$Gene,]

# sort the genes before transposing
# could probably do this more stably with melt + dcast...?
cat("sorting and transposing expressions...\n")
setorder(geuvadis.rnaseq.sub, Gene)
trnaseq = data.table(t(geuvadis.rnaseq.sub[,-c(1)]))
colnames(trnaseq) = geuvadis.rnaseq.sub$Gene
trnaseq = cbind(colnames(geuvadis.rnaseq.sub)[-1], trnaseq)
colnames(trnaseq)[1] = "SubjectID"

# put column names and a column of sample names
# make sure to sort by sample
setorder(trnaseq, SubjectID)

# melt the predicted and measured expression data.tables
cat("melting both prediction and measured expression data tables...\n")
pred.melted   = melt(geuvadis.predictions.sub, id.vars = c("Gene"))
rnaseq.melted = melt(geuvadis.rnaseq.sub, id.vars = "Gene")

# standardize their column names too
colnames(pred.melted)   = c("Gene", "SubjectID", "Predicted_Expr")
colnames(rnaseq.melted) = c("Gene", "SubjectID", "Measured_Expr")

# now merge the data frames
cat("merging melted data tables...\n")
geuvadis.rnapred = merge(pred.melted, rnaseq.melted, by = c("Gene", "SubjectID"))
setorder(geuvadis.rnapred, Gene, SubjectID)

# do two linear models
# first regress all (Gene, SubjectID) pairs of predicted expression onto measured expression
cat("first linear model: regress predicted onto measured expression for all (gene, subject) pairs\n") 
sink(out.lm.file)
print(summary(lm(Predicted_Expr ~ Measured_Expr, data = geuvadis.rnapred)))
sink()

# now do individual gene regressions
cat("second linear model: individual regressions for each gene\n")
compute.gene.lm(geuvadis.rnapred, out.genelm.file)

# now let us compare crosspopulation predictions 
cat("now analyzing crosspopulation predictions\n")
altpop.predictions = fread(altpop.pred.file)
altpop.expr        = fread(altpop.exprfile)

# melt the data.table for the measured expression data
altpop.expr.melt = melt(data = altpop.expr, id.vars = c("Gene"))
colnames(altpop.expr.melt) = c("Gene", "SubjectID", "Measured_Expr")

# rename column names for altpop.predictions
colnames(altpop.predictions) = c("Gene", "SubjectID", "Predicted_Expr")

# merge the measured and predicted expressions in the alternative population 
altpop.rnapred = merge(altpop.expr.melt, altpop.predictions, by = c("Gene", "SubjectID"))

# first model: regress measured, predicted expression in all (gene, subject) pairs
cat("first linear model in other population: regress predicted onto measured expression for all (gene, subject) pairs\n")  
sink(altpop.out.lm.file)
print(summary(lm(Predicted_Expr ~ Measured_Expr, data = altpop.rnapred)))
sink()

cat("second linear model in other population: individual regressions for each gene\n")
compute.gene.lm(altpop.rnapred, altpop.out.genelm.file)

# show any warnings
cat("any warnings?\n")
warnings()

# done!
cat("all done!\n\n")
