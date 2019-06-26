#!/usr/bin/env Rscript --vanilla

# ==========================================================================================
# load libraries
# ==========================================================================================
suppressMessages(library(methods))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(broom))
suppressMessages(library(optparse))

# ==========================================================================================
# subroutines
# ==========================================================================================
compute.r2.corr.onegene = function(x){
    the.fit = summary(lm(Predicted_Expr ~ Measured_Expr, data = x)) 
    my.cortest =  cor.test(~ Predicted_Expr + Measured_Expr, data = x, method = "spearman")
    my.output = data.frame(
        "Correlation" = my.cortest$estimate,
        "Corr_pval" = my.cortest$p.value,
        "R2" = the.fit$r.squared,
        "R2_pval" = the.fit$coefficients[2,4]
    )
    return(my.output)
}

compute.r2.corr = function(df){
    new.df = df %>% 
        na.omit %>% 
        dplyr::group_by(Train_Pop, Test_Pop, Gene) %>% 
        do(compute.r2.corr.onegene(.)) %>% 
        as.data.table
    return(new.df)
}

# ==========================================================================================
# parse command line arguments 
# ==========================================================================================

option_list = list(
    make_option(
        c("-a", "--EUR373-to-EUR373"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR373 to all Europeans",
        metavar = "character"
    ),
    make_option(
        c("-b", "--EUR373-to-AFR"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR373 to the Yoruba",
        metavar = "character"
    ),
    make_option(
        c("-c", "--EUR278-to-EUR278"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR278 to non-Finnish Europeans",
        metavar = "character"
    ),
    make_option(
        c("-d", "--EUR278-to-FIN"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR278 to Finns",
        metavar = "character"
    ),
    make_option(
        c("-e", "--EUR278-to-AFR"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR278 to the Yoruba",
        metavar = "character"
    ),
    make_option(
        c("-f", "--AFR-to-EUR373"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR278 to the Yoruba",
        metavar = "character"
    ),
    make_option(
        c("-g", "--AFR-to-AFR"),
        type    = "character",
        default = NULL,
        help    = "Predictions from AFR to the Yoruba",
        metavar = "character"
    ),
    make_option(
        c("-i", "--output-predictions"),
        type    = "character",
        default = NULL,
        help    = "Path to compiled data table of crosspopulation predictions",
        metavar = "character"
    ),
    make_option(
        c("-j", "--output-r2"),
        type    = "character",
        default = NULL,
        help    = "Path to compiled R2 results",
        metavar = "character"
    ),
    make_option(
        c("-k", "--EUR278-ids"),
        type    = "character",
        default = NULL,
        help    = "PLINK-formatted sample IDs for EUR278",
        metavar = "character"
    ),
    make_option(
        c("-l", "--FIN-ids"),
        type    = "character",
        default = NULL,
        help    = "PLINK-formatted sample IDs for FIN",
        metavar = "character"
    ),
    make_option(
        c("-o", "--poscorr"),
        type    = "character",
        default = NULL,
        help    = "Path to save list of genes in each train-test scenario with positive correlation between data and prediction",
        metavar = "character"
    ),
    make_option(
        c("-p", "--common-genes-poscorr"),
        type    = "character",
        default = NULL,
        help    = "Path to save list of genes with positive correlation between data and predictions in all train-test scenarios",
        metavar = "character"
    )
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n")
print(opt)


eur373.to.eur373.path = opt$EUR373_to_EUR373
eur373.to.afr.path    = opt$EUR373_to_AFR
eur278.to.eur278.path = opt$EUR278_to_EUR278
eur278.to.fin.path    = opt$EUR278_to_FIN
eur278.to.afr.path    = opt$EUR278_to_AFR
afr.to.eur373.path    = opt$AFR_to_EUR373
eur278.ids.path       = opt$EUR278_ids
fin.ids.path          = opt$FIN_ids
afr.to.afr.path       = opt$AFR_to_AFR
output.path.predictions = opt$output_predictions
output.path.r2        = opt$output_r2
output.path.r2.poscorr = opt$poscorr
output.path.r2.poscorr.commongenes = opt$common_genes_poscorr

# load data frames
eur373.to.eur373 = fread(eur373.to.eur373.path, header = TRUE)
eur373.to.afr    = fread(eur373.to.afr.path, header = TRUE)
eur278.to.eur278 = fread(eur278.to.eur278.path, header = TRUE)
eur278.to.fin    = fread(eur278.to.fin.path, header = TRUE)
eur278.to.afr    = fread(eur278.to.afr.path, header = TRUE)
afr.to.eur373    = fread(afr.to.eur373.path, header = TRUE)
afr.to.afr       = fread(afr.to.afr.path, header = TRUE)

eur278.ids = fread(eur278.ids.path, header = FALSE)
fin.ids    = fread(fin.ids.path , header = FALSE) 

# subset eur278, fin from AFR predictions
afr.to.eur278 = afr.to.eur373 %>% dplyr::filter(SubjectID %in% eur278.ids[[1]]) %>% as.data.table 
afr.to.fin    = afr.to.eur373 %>% dplyr::filter(SubjectID %in% fin.ids[[1]]) %>% as.data.table 

# make char vectors of train/test orders
train.pops = c("EUR373", "EUR373", "EUR278", "EUR278", "EUR278", "AFR", "AFR", "AFR", "AFR")
test.pops  = c("EUR373", "AFR", "EUR278", "FIN", "AFR", "EUR373", "EUR278", "FIN", "AFR")

# make list of data frames in same order
datatables = list(eur373.to.eur373, eur373.to.afr, eur278.to.eur278, eur278.to.fin, eur278.to.afr, afr.to.eur373, afr.to.eur278, afr.to.fin, afr.to.afr)

# add pop labels
n = length(train.pops)
for (i in 1:n) {
    datatables[[i]]$Train_Pop = train.pops[i]
    datatables[[i]]$Test_Pop = test.pops[i]
    datatables[[i]]$Spearman.rho = as.numeric(datatables[[i]]$Spearman.rho)
    setkey(datatables[[i]], Gene, P.value, R2, Spearman.rho, Train_Pop, Test_Pop)
}

# merge data frames together and save to file
geuvadis.predictions = Reduce(function(x,y) {merge(x, y, all = TRUE)}, datatables)
fwrite(x = geuvadis.predictions, file = output.path.predictions, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t") 

# summarize prediction results
geuvadis.results = compute.r2.corr(geuvadis.predictions)

# get genes with results in all populations
transcripts = geuvadis.results %>%
    group_by(Gene) %>%
    summarize(count = n()) %>%
    dplyr::filter(count == n) %>%
    select(Gene) %>%
    as.data.table

# use "transcripts" as a filtering variable to parse R2
r2 = geuvadis.results %>%
    dplyr::filter(Gene %in% transcripts[[1]]) %>%
    group_by(Train_Pop, Test_Pop) %>%
    summarize(r2 = mean(R2, na.rm = T), corr = mean(Spearman.rho, na.rm = T), Transcripts = n()) %>%
    as.data.table

# save filtered R2 to file
fwrite(x = r2, file = output.path.r2, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t") 

# create second filter against those transcripts, now only pulling positive correlations 
r2.poscorr = geuvadis.results %>%
    dplyr::filter(Gene %in% transcripts[[1]] & Spearman.rho > 0) %>%
    group_by(Train_Pop, Test_Pop) %>%
    summarize(r2 = mean(R2, na.rm = TRUE), corr = mean(Spearman.rho, na.rm = TRUE), Transcripts = n()) %>%
    as.data.table

# save refiltered r2 to file
fwrite(x = r2.poscorr, file = output.path.r2.poscorr, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t") 

# get genes that have a positive correlation in EACH population
transcripts.poscorr = geuvadis.results %>%
    dplyr::filter(Spearman.rho > 0) %>%
    group_by(Gene) %>%
    summarize(count = n()) %>%
    dplyr::filter(count == n) %>%
    select(Gene) %>%
    as.data.table

r2.poscorr.commongenes = geuvadis.results %>%
    dplyr::filter(Gene %in% transcripts.poscorr[[1]]) %>%
    group_by(Train_Pop, Test_Pop) %>%
    summarize(r2 = mean(R2, na.rm = TRUE), corr = mean(Spearman.rho, na.rm = TRUE), Transcripts = n()) %>%
    as.data.table

# save this last list to file
fwrite(x = r2.poscorr.commongenes, file = output.path.r2.poscorr.commongenes, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t") 

# want to know if distributions between poscorr genes, all genes are statistically significant
# must group by poscorr status (0: not poscorr gene, 1: poscorr gene) and continental origin (EUR or AFR)
poscorr.dunn.test.results = geuvadis.results %>%
    mutate(
        Train_Test = paste(Train_Pop, Test_Pop, sep = "_"),
        poscorr = as.integer(Gene %in% transcripts.poscorr[[1]]),
        Train_Continent = ifelse(Train_Pop == "AFR", "AFR", "EUR"),
        Test_Continent = ifelse(Test_Pop == "AFR", "AFR", "EUR"),
        Train_Test_Continent = paste(Train_Continent, Test_Continent, sep = "_")
    ) %>%
    dunn.test(
        x = with(., R2),
        g = with(., interaction(poscorr, Train_Test_Continent)),
        method = "bonferroni",
        kw = TRUE,
        label = TRUE
    ) %>% as.data.frame %>% as.data.table

# save these results to file
# for manuscript, the following comparisons are particularly relevant:
# 0.EUR_EUR - 1.EUR_EUR (predictions within EUR, comparing in/not in poscorr genes)
# 0.AFR_AFR - 1.AFR_AFR (predictions within AFR, *...*)
fwrite(x = poscorr.dunn.test.results, file = paste(output.path.r2.poscorr.commongenes, "dunntest.results", sep = "."), sep = "\t", quote = FALSE)
