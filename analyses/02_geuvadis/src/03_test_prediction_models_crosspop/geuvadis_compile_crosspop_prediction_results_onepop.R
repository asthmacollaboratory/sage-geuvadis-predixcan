library(optparse)
options(warn = -1)

# ========================================================================================
# parse command line variables
# ========================================================================================
option_list = list(
    make_option(
        c("-p", "--pop"),
        type    = "character",
        default = NULL, 
        help    = "Name of the population used for building prediction models",
        metavar = "character"
    ),
    make_option(
        c("-a", "--altpop"),
        type    = "character",
        default = NULL, 
        help    = "Name of the population used for validating prediction models",
        metavar = "character"
    ),
    make_option(
        c("-rp", "--rna-pop"),
        type    = "character",
        default = NULL, 
        help    = "RNA data from the population used for building prediction models",
        metavar = "character"
    ),
    make_option(
        c("-ra", "--rna-altpop"),
        type    = "character",
        default = NULL, 
        help    = "RNA data from the population used for validating prediction models",
        metavar = "character"
    ),
    make_option(
        c("-pp", "--predictions-pop"),
        type    = "character",
        default = NULL, 
        help    = "Predictions from the population used for building prediction models",
        metavar = "character"
    ),
    make_option(
        c("-pa", "--predictions-altpop"),
        type    = "character",
        default = NULL, 
        help    = "Predictions from the population used for validating prediction models",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-prefix"),
        type    = "character",
        default = NULL, 
        help    = "Prefix for output files",
        metavar = "character"
    ),
    make_option(
        c("-k", "--sample-pop-key"),
        type    = "character",
        default = NULL, 
        help    = "Path to two-column tab-separated file with SAMPLE and POP code",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

pop    = opt$pop
altpop = opt$pop

pop.rna.path  = opt$rna_pop
pop.pred.path = opt$predictions_pop

altpop.rna.path  = opt$rna_altpop
altpop.pred.path = opt$predictions_altpop

output.pfx = opt$output_prefix

key.file.path = opt$sample_pop_key

# =======================================================================================
# load libraries
# =======================================================================================
library(methods)
library(data.table)
library(dplyr)
library(purrr)
library(broom)
library(stringr)

# =======================================================================================
# subroutines
# =======================================================================================

# subroutine for extracting both Spearman rho, p-value from correlation test
cortest = function(x) {
    my.cortest =  cor.test(~ Predicted_Expr + Measured_Expr, data = x, method = "spearman")
    return(data.table(Gene = x$Gene, Correlation = my.cortest$estimate, Corr.p.value = my.cortest$p.value))
}

compute.r2.corr.onegene = function(x){
    the.fit = summary(lm(Predicted_Expr ~ Measured_Expr, data = x)) 
    my.cortest =  cor.test(~ Predicted_Expr + Measured_Expr, data = x, method = "spearman")
    my.output = data.frame("Correlation" = my.cortest$estimate, "Corr_pval" = my.cortest$p.value, "R2" = the.fit$r.squared, "R2_pval" = the.fit$coefficients[2,4])
    return(my.output)
}

compute.r2.corr = function(df){
    new.df = df %>% 
        na.omit %>% 
        dplyr::group_by(Gene, Test_Pop) %>% 
        do(compute.r2.corr.onegene(.)) %>% 
        as.data.table
    return(new.df)
}


# ========================================================================================
# load files 
# ========================================================================================

pop.rna  = fread(pop.rna.path)
pop.pred = fread(pop.pred.path)

altpop.rna  = fread(altpop.rna.path)
altpop.pred = fread(altpop.pred.path)

# key.file structure; no header, 1st col Subject ID, 2nd Pop
key.file = fread(key.file.path, header = FALSE)
colnames(key.file) = c("SubjectID", "Test_Pop") 

# melt rna files since they are arranged as matrices
pop.rna.melt    = melt(pop.rna, id.vars = "Gene", variable.name = "SubjectID")
altpop.rna.melt = melt(altpop.rna, id.vars = "Gene", variable.name = "SubjectID")
colnames(pop.rna.melt)[3]    = "Measured_Expr"
colnames(altpop.rna.melt)[3] = "Measured_Expr"

# melt predictions if necessary
if (dim(pop.pred)[2] != 3){
    pop.pred = melt(pop.pred, id.vars = "Gene", variable.name = "SubjectID")
}
colnames(pop.pred)[3] = "Predicted_Expr"

if (dim(altpop.pred)[2] != 3){
    altpop.pred = melt(altpop.pred, id.vars = "Gene", variable.name = "SubjectID")
}
colnames(altpop.pred)[3] = "Predicted_Expr"

# merge predictions and measurements
pop.rnapred    = merge(pop.rna.melt, pop.pred, by = c("Gene", "SubjectID"))
altpop.rnapred = merge(altpop.rna.melt, altpop.pred, by = c("Gene", "SubjectID"))

# merge rnapred data.tables
# we will keep populations separate using key.file later
rnapred = rbind(pop.rnapred, altpop.rnapred)

# add population code
# use a right_join to ensure that we only save predictions in crosspop scheme
# this discards the samples not in the 89 selected for each pop
# result is a data.frame, so recast it to data.table
rnapred = as.data.table(right_join(rnapred, key.file, by = "SubjectID"))
rnapred$Train_Pop = str_to_upper(pop) 

# recover memory
pop.rna         = FALSE
pop.pred        = FALSE
altpop.rna      = FALSE
altpop.pred     = FALSE
pop.rna.melt    = FALSE
altpop.rna.melt = FALSE
pop.pred        = FALSE
altpop.pred     = FALSE
pop.rnapred     = FALSE
altpop.rnapred  = FALSE
pop.rna.melt    = FALSE
gc()


allpop.results = compute.r2.corr(rnapred) 
allpop.results$Train_Pop = str_to_upper(pop)

# save results
allpop.results.path = paste0(output.pfx, ".", pop, ".predictinto.allpop.results.txt")
fwrite(x = allpop.results, file = allpop.results.path, quote = FALSE, na = "NA", sep = "\t")

# want to compute pop-wise summaries of these results
# first get the predicted genes in common across all pops
genes.in.common = allpop.results %>%
    na.omit %>%
    count(Gene) %>%
    dplyr::filter(., n == 5) %>% ### TODO: determine number of pops programmatically? currently set to 5
    distinct %>%
    select(Gene) %>%
    as.data.table

# save these genes to file
commongenes.path = paste0(output.pfx, ".", pop, ".commongenes.txt")
fwrite(x = genes.in.common, file = commongenes.path, quote = FALSE, sep = "\t", col.names = FALSE)

# get mean R2
# first computes mean R2 for each gene in each pop
# then reduces across genes to get average R2 for pop
mean.r2 = allpop.results %>%
    select(Gene, Train_Pop, Test_Pop, R2) %>%
    dplyr::filter(., Gene %in% genes.in.common$Gene) %>%
    group_by(Gene, Test_Pop) %>%
    summarize(R2_mean = mean(R2, na.rm = TRUE)) %>%
    group_by(Test_Pop) %>%
    summarize(Mean = mean(R2_mean, na.rm = TRUE),
              Std_Err = sd(R2_mean, na.rm = TRUE),
              Num_Genes = n()
    ) %>%
    as.data.table

# rename columns to clarify the training population
colnames(mean.r2)[2] = "R2_Mean" 
colnames(mean.r2)[3] = "R2_SD"
mean.r2$Train_Pop = str_to_upper(pop)

# reorder columns
mean.r2 = mean.r2 %>% select(Train_Pop, Test_Pop, R2_Mean, R2_SD, Num_Genes)

# save to file
allpop.summary.path = paste0(output.pfx, ".", pop, ".predictinto.allpop.summary.txt")
fwrite(x = mean.r2, file = allpop.summary.path, quote = FALSE, na = "NA", sep = "\t")

# repeat search of genes in common across pops
# this time purge those with nonpositive correlation
genes.in.common.poscorr = allpop.results %>%
    dplyr::filter(., Correlation > 0) %>% # <-- operative difference
    na.omit %>%
    count(Gene) %>%
    dplyr::filter(., n == 5) %>%
    distinct %>%
    select(Gene) %>%
    as.data.table

# save these genes to file
commongenes.path.poscorr = paste0(output.pfx, ".", pop, ".commongenes.poscorr.txt")
fwrite(x = genes.in.common.poscorr, file = commongenes.path.poscorr, quote = FALSE, sep = "\t", col.names = FALSE)

# get another mean R2
# this one discards negative correlations
# produces more realistic R2s but discards a lot of genes
mean.r2.poscorr = allpop.results %>%
    dplyr::filter(., Correlation > 0) %>%
    select(Gene, Test_Pop, R2) %>%
    group_by(Gene) %>%
    dplyr::filter(., Gene %in% genes.in.common.poscorr$Gene) %>%
    group_by(Gene, Test_Pop) %>%
    summarize(R2_mean = mean(R2, na.rm = TRUE)) %>%
    group_by(Test_Pop) %>%
    summarize(Mean = mean(R2_mean, na.rm = TRUE),
              Std_Err = sd(R2_mean, na.rm = TRUE),
              Num_Genes = n()
    ) %>%
    as.data.table

# rename columns to clarify the training population
colnames(mean.r2.poscorr)[2] = "R2_Mean" 
colnames(mean.r2.poscorr)[3] = "R2_SD"
mean.r2.poscorr$Train_Pop = str_to_upper(pop)

# reorder columns
mean.r2.poscorr = mean.r2.poscorr %>% select(Train_Pop, Test_Pop, R2_Mean, R2_SD, Num_Genes)

# save to file
allpop.summary.path.poscorr = paste0(output.pfx, ".", pop, ".predictinto.allpop.summary.poscorr.txt")
fwrite(x = mean.r2.poscorr, file = allpop.summary.path.poscorr, quote = FALSE, na = "NA", sep = "\t")
