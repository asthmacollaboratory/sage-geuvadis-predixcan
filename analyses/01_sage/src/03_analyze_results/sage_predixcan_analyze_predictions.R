#!/usr/bin/env Rscript --vanilla
# =======================================================================================
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys
#
# This script processes and analyzes PrediXcan predictions in SAGE data 
# =======================================================================================

# =======================================================================================
# load libraries 
# =======================================================================================
library(data.table)
library(purrr)
library(broom)
library(ggplot2)
library(dplyr)
library(optparse)

# parse command line arguments
option_list = list(
    make_option(
        c("-a", "--predixcan-dir"),
        type    = "character",
        default = NULL,
        help    = "Data where PrediXcan output are stored.",
        metavar = "character"
    ),
    make_option(
        c("-b", "--expression-file"),
        type    = "character",
        default = NULL,
        help    = "File where RNA-Seq data are stored.",
        metavar = "character"
    ),
    make_option(
        c("-c", "--output-dir"),
        type    = "character",
        default = NULL,
        help    = "Directory to store output.",
        metavar = "character"
    ),
    make_option(
        c("-d", "--plot-type"),
        type    = "character",
        default = "png",
        help    = "Desired file type for plots of results, given as a file extension (e.g. PNG, JPEG, PDF, TIFF, SVG)",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

# =======================================================================================
# subroutines 
# =======================================================================================

lmtest = function(x, lm.formula) {
    x.nona = na.omit(x)
    n = dim(x.nona)[1]
    if ( n < 1 ) {
        return(data.table(Gene = x$Gene, R2 = NA, N = 0))
    } 
    my.lm = lm(formula(lm.formula), data = x.nona, na.action = na.omit)
    return(data.table(Gene = x$Gene, R2 = summary(my.lm)$r.squared, N = n))
}

# subroutine for extracting both Spearman rho, p-value from correlation test
cortest = function(x, cor.formula) {
    x.nona = na.omit(x)
    n = dim(x.nona)[1]
    if ( n < 1 ) {
        return(data.table(Gene = x$Gene, Correlation = NA, Corr.p.value = NA))
    } 
    my.cortest = cor.test(formula(cor.formula), x.nona, method = "spearman", na.action = na.omit)
    return(data.table(Gene = x.nona$Gene, Correlation = my.cortest$estimate, Corr.p.value = my.cortest$p.value))
}

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
        dplyr::group_by(Prediction.Weights, Gene) %>% 
        do(compute.r2.corr.onegene(.)) %>% 
        as.data.table
    return(new.df)
}


plot.boxplot.r2 = function(x, ngenes, scale.fill.palette, boxplot.labels) {

    my.boxplot = ggplot(x, aes(x = Prediction.Weights, y = R2, fill = Prediction.Weights)) +
        geom_boxplot() +
        xlab("Prediction Weight Set") +
        ylab(expression(R^{2})) +
        ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ngenes) ~ "genes")) +
        scale_fill_manual(values = scale.fill.palette) +
        scale_x_discrete(labels = boxplot.labels) +
        ylim(-0.1, 1) +
        theme(
            legend.position = "none",
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12)
        )

    return(my.boxplot)
}


plot.violinplot.r2 = function(x, ngenes) {

    my.violinplot = ggplot(x, aes(x = Prediction.Weights, y = R2)) +
        geom_violin() +
        xlab("Prediction Weight Set") +
        ylab(expression(R^{2})) +
        ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ngenes) ~ "genes"))

    return(my.violinplot)
}

plot.hist.r2 = function(x, ngenes) {

    my.hist = ggplot(x, aes(x = R2)) +
        geom_histogram(aes(y = ..density..)) +
        geom_density() +
        xlab("Prediction Weight Set") +
        ylab(expression(R^{2})) +
        ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ ngenes ~ "genes")) +
        facet_wrap(~ Prediction.Weights)
 
    return(my.hist)
}


save.plots = function(g1, name1, g2, name2, g3, name3){
    # boxplot
    ggsave(plot = g1, filename = name1, dpi = 300, type = "cairo", width = 15, height = 5, units = "in")

    # histogram
    ggsave(plot = g2, filename = name2, dpi = 300, type = "cairo")

    # violinplot
    ggsave(plot = g3, filename = name3, dpi = 300, type = "cairo")

    return()
}

# =======================================================================================
# file and directory paths 
# =======================================================================================
predixcan.dir  = opt$predixcan_dir 
#data.dir       = opt$data.dir
output.dir     = opt$output_dir
plot.type      = opt$plot_type

dgn.path       = file.path(predixcan.dir, "DGN",  "DGN_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
gtex6.path     = file.path(predixcan.dir, "GTEx", "GTEx_v6p_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
gtex7.path     = file.path(predixcan.dir, "GTEx", "GTEx_v7_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
mesa.afa.path  = file.path(predixcan.dir, "MESA", "MESA_AFA_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
mesa.afhi.path = file.path(predixcan.dir, "MESA", "MESA_AFHI_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
mesa.cau.path  = file.path(predixcan.dir, "MESA", "MESA_CAU_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")
mesa.all.path  = file.path(predixcan.dir, "MESA", "MESA_ALL_predixcan_prediction_ALLCHR_predicted_expression_melted.txt")

#sage.rna.path  = file.path(data.dir, "sage_39_wgs_for_rnaseq_expression_melted.txt")
sage.rna.path  = opt$expression_file

rdata.path     = file.path(output.dir, "sage_predixcan_allresults_allplots.Rdata")

# =======================================================================================
# load data
# =======================================================================================
dgn       = fread(dgn.path)
gtex6     = fread(gtex6.path)
gtex7     = fread(gtex7.path)
mesa.afa  = fread(mesa.afa.path)
mesa.afhi = fread(mesa.afhi.path)
mesa.cau  = fread(mesa.cau.path)
mesa.all  = fread(mesa.all.path)

sage = fread(sage.rna.path, header = TRUE)

# =======================================================================================
# merge data into single data frame  
# =======================================================================================

# must rename columns of each data.table, particularly the ones with predicted expression values
# this facilitates merging them later
repos = c("DGN", "GTEx_v6p", "GTEx_v7", "MESA_AFA", "MESA_AFHI", "MESA_CAU", "MESA_ALL")
repo.results = list(dgn, gtex6, gtex7, mesa.afa, mesa.afhi, mesa.cau, mesa.all)
for (i in 1:length(repos)) {
    colnames(repo.results[[i]]) = c("SubjectID", "Gene", "Predicted_Expr")
    repo.results[[i]]$Prediction.Weights =  repos[i]
}

# perform a full join of all prediction results
predixcan.all = repo.results %>% reduce(full_join, by = c("SubjectID", "Gene", "Prediction.Weights", "Predicted_Expr")) %>% as.data.table

# if necessary, melt measurements prior to merge
# when data.table 'sage. doesn't have 3 cols, we assume that it is in a BED format
if (ncol(sage) != 3) {
    sage.melt = melt(sage[,-c(1:3)], id.vars = "Gene", variable.name = "SubjectID", value.name = "Measured_Expr")
} else {
    sage.melt = sage
    names(sage.melt) = c("Gene", "SubjectID", "Measured_Expr")
}

# merge predictions with measurements
sage.predixcan.all = merge(sage.melt, predixcan.all, by = c("Gene", "SubjectID"), all = TRUE)

# recover some memory
predixcan.all = sage = sage.melt = FALSE
gc()


# =======================================================================================
# analyze predictions 
# =======================================================================================

sage.predixcan.all.results = compute.r2.corr(sage.predixcan.all)
sage.predixcan.all = FALSE
gc()


# =======================================================================================
# plot results
# =======================================================================================

# gplot2 colorblind-friendly palette with grey:
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73") 

# set boxplot labels
my.boxplot.labels = c("GTEx v6p", "GTEx v7", "DGN", "MESA_AFA", "MESA_AFHI", "MESA_CAU", "MESA_ALL", "Test: GTEx v7", "Test: MESA AFA", "Test: MESA AFHI", "Test: MESA CAU", "Test: MESA ALL")

# compare imputation performance from all four prediction weights
# must first load GTEx v7 testing R2s from PredictDB
predixcan.gtex7.metrics.path     = file.path(predixcan.dir, "gtex7.test.metrics.txt")
predixcan.mesa.afa.metrics.path  = file.path(predixcan.dir, "mesa.AFA.test.metrics.txt")
predixcan.mesa.afhi.metrics.path = file.path(predixcan.dir, "mesa.AFHI.test.metrics.txt")
predixcan.mesa.cau.metrics.path  = file.path(predixcan.dir, "mesa.CAU.test.metrics.txt")
predixcan.mesa.all.metrics.path  = file.path(predixcan.dir, "mesa.ALL.test.metrics.txt")

predixcan.gtex7.metrics     = fread(predixcan.gtex7.metrics.path, header = TRUE) 
predixcan.mesa.afa.metrics  = fread(predixcan.mesa.afa.metrics.path, header = TRUE) 
predixcan.mesa.afhi.metrics = fread(predixcan.mesa.afhi.metrics.path, header = TRUE) 
predixcan.mesa.cau.metrics  = fread(predixcan.mesa.cau.metrics.path, header = TRUE) 
predixcan.mesa.all.metrics  = fread(predixcan.mesa.all.metrics.path, header = TRUE) 

# must rename columns of each data.table, particularly the ones with predicted expression values
# this facilitates merging them later
repos.test = c("Test: GTEx_v7", "Test: MESA_AFA", "Test: MESA_AFHI", "Test: MESA_CAU", "Test: MESA_ALL")
repo.test = list(predixcan.gtex7.metrics, predixcan.mesa.afa.metrics, predixcan.mesa.afhi.metrics, predixcan.mesa.cau.metrics, predixcan.mesa.all.metrics)
for (i in 1:length(repos.test)) {

    # note that at this point the 2nd column actually contains HUGO identifiers for genes
    colnames(repo.test[[i]]) = c("Gene", "Prediction.Weights", "R2", "Correlation")

    # 2nd column now overwritten with repo name
    repo.test[[i]]$Prediction.Weights = repos.test[i]

    # trim .X, the transcript number to the ENSG ID, starting at 15th character
    repo.test[[i]]$Gene = strtrim(repo.test[[i]]$Gene, 15)
}

# perform a full join of all prediction results
predixcan.all.test = rbindlist(repo.test)
repo.test = FALSE; gc()  ## recover memory

# extract coefficients of determination (R2) from genewise summary results
all.results = sage.predixcan.all.results %>%
    select(Prediction.Weights, Gene, R2, Correlation) %>%
    as.data.table %>%
    rbind(., predixcan.all.test, use.names = TRUE, fill = TRUE) 

r2.all = all.results %>%
    select(Gene, Prediction.Weights, R2) %>%
    as.data.table %>%
    na.omit
colnames(r2.all) = c("Gene", "Prediction.Weights", "R2")
setorderv(r2.all, "Prediction.Weights")
ngenes.r2 = r2.all %>% select(Gene) %>% unlist %>% sort %>% unique %>% length

# will produce three kinds of summary plots:
# (1): boxplot
# (2): violin plot
# (3): histogram + density plot 
my.boxplot.r2 = ggplot(r2.all, aes(x = Prediction.Weights, y = R2, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ngenes.r2) ~ "genes")) +
    scale_fill_manual(values = cbPalette) +
    scale_x_discrete(labels = my.boxplot.labels) +
    ylim(-0.1, 1) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
    )

my.violinplot.r2 = ggplot(r2.all, aes(x = Prediction.Weights, y = R2)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ngenes.r2) ~ "genes"))
my.hist.r2 = ggplot(r2.all, aes(x = R2)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ngenes.r2) ~ "genes")) +
    facet_wrap(~ Prediction.Weights)

#my.boxplot.r2 = plot.boxplot.r2(r2.all, ngenes.r2, cbPalette, my.boxplot.labels)
#my.violinplot.r2 = plot.violinplot.r2(r2.all, ngenes.r2)
#my.hist.r2 = plot.hist.r2(r2.all, ngenes.r2)

# save plots to file
boxplot.r2.path = file.path(output.dir, paste0("sage.predixcan.r2.all.boxplot.", plot.type))
violinplot.r2.path = file.path(output.dir, paste0("sage.predixcan.r2.all.violinplot.", plot.type))
histogram.r2.path = file.path(output.dir, paste0("sage.predixcan.r2.all.histogram.", plot.type))

save.plots(my.boxplot.r2, boxplot.r2.path, my.violinplot.r2, violinplot.r2.path, my.hist.r2, histogram.r2.path)


# do for only common genes
commongenes.r2 = r2.all %>%
    count(Gene) %>%
    dplyr::filter(n > 11) %>%
    select(Gene) %>%
    as.data.table
r2.commongenes = r2.all %>% dplyr::filter(Gene %in% commongenes.r2$Gene)
ncommongenes.r2 = r2.commongenes %>% select(Gene) %>% unlist %>% sort %>% unique %>% length 

my.boxplot.r2.commongenes = ggplot(r2.commongenes, aes(x = Prediction.Weights, y = R2, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of " ~ R^2 ~ "across different prediction weight sets over" ~ .(ncommongenes.r2) ~ "common genes")) +
    scale_fill_manual(values = cbPalette) +
    scale_x_discrete(labels = my.boxplot.labels) +
    ylim(-0.1, 1) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
    )

my.violinplot.r2.commongenes = ggplot(r2.commongenes, aes(x = Prediction.Weights, y = R2)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ncommongenes.r2) ~ "common genes"))

my.hist.r2.commongenes = ggplot(r2.commongenes, aes(x = R2)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ncommongenes.r2) ~ "common genes")) +
    facet_wrap(~ Prediction.Weights)

#my.boxplot.r2.commongenes = plot.boxplot.r2(commongenes.r2, ncommongenes.r2, cbPalette, my.boxplot.labels)
#my.violinplot.r2.commongenes = plot.violinplot.r2(commongenes.r2, ncommongenes.r2)
#my.hist.r2.commongenes = plot.hist.r2(commongenes.r2, ncommongenes.r2)

# save plots to file
boxplot.r2.commongenes.path = file.path(output.dir, paste0("sage.predixcan.r2.commongenes.boxplot.", plot.type))
violinplot.r2.commongenes.path = file.path(output.dir, paste0("sage.predixcan.r2.commongenes.violinplot.", plot.type))
histogram.r2.commongenes.path = file.path(output.dir, paste0("sage.predixcan.r2.commongenes.histogram.", plot.type))

save.plots(my.boxplot.r2.commongenes, boxplot.r2.commongenes.path, my.violinplot.r2.commongenes, violinplot.r2.commongenes.path, my.hist.r2.commongenes, histogram.r2.commongenes.path)


# do same, but for correlations instead of R2s
rho.all = all.results %>% 
    select(Gene, Prediction.Weights, Correlation) %>% 
    as.data.table %>%
    na.omit
#colnames(rho.all) = c("Gene", "Prediction.Weights", "Correlation")
setorderv(rho.all, "Prediction.Weights")
ngenes.rho = ngenes.r2

my.hist.rho = ggplot(rho.all, aes(x = Correlation)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab("Spearman rho") +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ngenes.rho, " genes")) +
    facet_wrap(~ Prediction.Weights)
my.boxplot.rho = ggplot(rho.all, aes(x = Prediction.Weights, y = Correlation, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(paste("Correlations (Spearman ", rho, ")"))) +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ngenes.rho, " genes")) +
    scale_fill_manual(values = cbPalette) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
    ) +
    scale_x_discrete(labels = my.boxplot.labels)
my.violinplot.rho = ggplot(rho.all, aes(x = Prediction.Weights, y = Correlation)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab("Correlations") +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ngenes.rho, " genes"))

boxplot.rho.path = file.path(output.dir, paste0("sage.predixcan.rho.all.boxplot.", plot.type))
violinplot.rho.path = file.path(output.dir, paste0("sage.predixcan.rho.all.violinplot.", plot.type))
histogram.rho.path = file.path(output.dir, paste0("sage.predixcan.rho.all.histogram.", plot.type))

save.plots(my.boxplot.rho, boxplot.rho.path, my.violinplot.rho, violinplot.rho.path, my.hist.rho, histogram.rho.path)

# could filter for rho in common
# but for basis of comparison, we really should focus on common genes from R2
#rho.commongenes = all.rho %>%
#    group_by(Gene) %>%
#    tally %>%
#    dplyr::filter(n > 11) %>%
#    select(Gene) %>%
#    as.data.table
#all.rho.commongenes = all.rho %>% dplyr::filter(Gene %in% rho.commongenes$Gene)
#ncommongenes.rho = all.rho.commongenes %>% select(Gene) %>% unlist %>% sort %>% unique %>% length 
rho.commongenes  = rho.all %>% filter(Gene %in% commongenes.r2$Gene)
commongenes.rho  = commongenes.r2
ncommongenes.rho = ncommongenes.r2


my.hist.rho.commongenes = ggplot(rho.commongenes, aes(x = Correlation)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab("Spearman rho") +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ncommongenes.rho, " genes in common")) +
    facet_wrap(~ Prediction.Weights)
my.boxplot.rho.commongenes = ggplot(rho.commongenes, aes(x = Prediction.Weights, y = Correlation, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(paste("Correlations (Spearman ", rho, ")"))) +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ncommongenes.rho, " genes in common")) +
    scale_fill_manual(values = cbPalette) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
    ) +
    scale_x_discrete(labels = my.boxplot.labels)
my.violinplot.rho.commongenes = ggplot(rho.commongenes, aes(x = Prediction.Weights, y = Correlation)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab("Correlations") +
    ggtitle(paste0("Distribution of correlations across different prediction weight sets for ", ncommongenes.rho, " genes in common"))

boxplot.rho.commongenes.path = file.path(output.dir, paste0("sage.predixcan.rho.commongenes.boxplot.", plot.type))
violinplot.rho.commongenes.path = file.path(output.dir, paste0("sage.predixcan.rho.commongenes.violinplot.", plot.type))
histogram.rho.commongenes.path = file.path(output.dir, paste0("sage.predixcan.rho.commongenes.histogram.", plot.type))

save.plots(my.boxplot.rho.commongenes, boxplot.rho.commongenes.path, my.violinplot.rho.commongenes, violinplot.rho.commongenes.path, my.hist.rho.commongenes, histogram.rho.commongenes.path)


# can now compile some summaries

# how many genes predicted + measured in each repo?
genes.predicted.perrepo = all.results %>% na.omit %>% count(Prediction.Weights) %>% as.data.table

# how many genes predicted + measured with positive correlation in each repo?
genes.predicted.poscorr.perrepo = all.results %>% na.omit %>% filter(Correlation > 0) %>% count(Prediction.Weights) %>% as.data.table

# average r2 by repo?
r2.summaries  = r2.all %>% group_by(Prediction.Weights) %>% summarize(r2 = mean(R2, na.rm = T)) %>% as.data.table

# average r2 by repo over common genes?
r2.commongenes.summaries  = r2.commongenes %>% group_by(Prediction.Weights) %>% summarize(r2 = mean(R2, na.rm = T)) %>% as.data.table

# average rho by repo?
rho.summaries = rho.all %>% group_by(Prediction.Weights) %>% summarize(rho = mean(Correlation, na.rm = T)) %>% as.data.table

# average rho by repo over common genes?
rho.commongenes.summaries = rho.commongenes %>% group_by(Prediction.Weights) %>% summarize(rho = mean(Correlation, na.rm = T)) %>% as.data.table

# last detail:
# make plots for common genes with poscorr
# must get correlations for this to work
#all.r2.corr = merge(all.r2, all.rho, by = c("Gene", "Prediction.Weights"), all = TRUE)
commongenes.poscorr.r2 = all.results %>% 
    na.omit %>% 
    dplyr::filter(Correlation > 0) %>% 
    count(Gene) %>% 
    dplyr::filter(n > 11) %>% 
    select(Gene) %>% 
    as.data.table
#colnames(commongenes.r2) = c("Gene", "Prediction.Weights", "R2")
r2.commongenes.poscorr = r2.all %>% dplyr::filter(Gene %in% commongenes.poscorr.r2$Gene)
ncommongenes.r2.poscorr = r2.commongenes.poscorr %>% select(Gene) %>% unlist %>% unname %>% sort %>% unique %>% length 

my.boxplot.r2.commongenes = ggplot(r2.commongenes.poscorr, aes(x = Prediction.Weights, y = R2, fill = Prediction.Weights)) +
    geom_boxplot() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of " ~ R^2 ~ "across different prediction weight sets over" ~ .(ncommongenes.r2.poscorr) ~ "common genes")) +
    scale_fill_manual(values = cbPalette) +
    scale_x_discrete(labels = my.boxplot.labels) +
    ylim(-0.1, 1) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
    )

my.violinplot.r2.commongenes = ggplot(r2.commongenes.poscorr, aes(x = Prediction.Weights, y = R2)) +
    geom_violin() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ncommongenes.r2.poscorr) ~ "common genes"))

my.hist.r2.commongenes = ggplot(r2.commongenes.poscorr, aes(x = R2)) +
    geom_histogram(aes(y = ..density..)) +
    geom_density() +
    xlab("Prediction Weight Set") +
    ylab(expression(R^{2})) +
    ggtitle(bquote("Distribution of" ~ R^2 ~ "across different prediction weight sets over" ~ .(ncommongenes.r2.poscorr) ~ "common genes")) +
    facet_wrap(~ Prediction.Weights)

# save plots to file
boxplot.r2.commongenes.path = file.path(output.dir, paste0("sage.predixcan.r2.commongenes.poscorr.boxplot.", plot.type))
violinplot.r2.commongenes.path = file.path(output.dir, paste0("sage.predixcan.r2.commongenes.poscorr.violinplot.", plot.type))
histogram.r2.commongenes.path = file.path(output.dir, paste0("sage.predixcan.r2.commongenes.poscorr.histogram.", plot.type))

save.plots(my.boxplot.r2.commongenes, boxplot.r2.commongenes.path , my.violinplot.r2.commongenes, violinplot.r2.commongenes.path, my.hist.r2.commongenes, histogram.r2.commongenes.path)

# last plot: compare R2 in SAGE and GTEx 7
# will filter data into format easy for plotting
# will throw in summary since those numbers appear in text
sage.gtex7.r2 = r2.all %>%
    # extract just GTEx7 info
    filter(Prediction.Weights %in% c("GTEx_v7", "Test: GTEx_v7")) %>% 
    # recast into N x 2 data.frame for easier plotting
    dcast(., Gene ~ Prediction.Weights) %>%
    # focus on GTEx v7 genes with testing R2 > 0.2
    filter(`Test: GTEx_v7` > 0.2) %>%
    as.data.table %>%
    # rename the test R2 column for easier plotting
    rename(., "Test_R2" = "Test: GTEx_v7")

# summary
summary(sage.gtex7.r2)

# now construct plot
sage.gtex7.plot = sage.gtex7.r2 %>%
    ggplot(., aes(x = Test_R2, y = GTEx_v7)) +
        geom_point(color = "black") +
        # add regression line (blue) to compare points
        geom_smooth(se = TRUE, method = "lm") +
        # add y=x line (red dotted)
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        ggtitle("PredictDB performance in SAGE with GTEx v7 weights") +
        xlab(bquote("PredictDB training "~ R^2 ~ " in GTEx v7")) +
        ylab(bquote("PrediXcan " ~ R^2 ~ " in SAGE"))


# repeat for MESA
sage.mesa.all.r2 = r2.all %>%
    filter(Prediction.Weights %in% c("MESA_ALL", "Test: MESA_ALL")) %>%
    dcast(., Gene ~ Prediction.Weights) %>%
    filter(`Test: MESA_ALL` > 0.2) %>%
    as.data.table %>%
    rename(., "Test_R2" = "Test: MESA_ALL")

sage.mesa.all.plot = sage.mesa.all.r2 %>% 
    ggplot(., aes(x = Test_R2, y = MESA_ALL)) +
        geom_point(color = "darkred") +
        # add regression line (blue) to compare points
        geom_smooth(se = TRUE, method = "lm") +
        # add y=x line (red dotted)
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        ggtitle("PredictDB performance in SAGE with MESA_ALL weights") +
        xlab(bquote("PredictDB training "~ R^2 ~ " in MESA_ALL")) +
        ylab(bquote("PrediXcan " ~ R^2 ~ " in SAGE"))

sage.mesa.afa.r2 = r2.all %>%
    filter(Prediction.Weights %in% c("MESA_AFA", "Test: MESA_AFA")) %>%
    dcast(., Gene ~ Prediction.Weights) %>%
    filter(`Test: MESA_AFA` > 0.2) %>%
    as.data.table %>%
    rename(., "Test_R2" = "Test: MESA_AFA")

sage.mesa.afa.plot = sage.mesa.afa.r2 %>% 
    ggplot(., aes(x = Test_R2, y = MESA_AFA)) +
        geom_point(color = "darkgreen") +
        # add regression line (blue) to compare points
        geom_smooth(se = TRUE, method = "lm") +
        # add y=x line (red dotted)
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        ggtitle("PredictDB performance in SAGE with MESA_AFA weights") +
        xlab(bquote("PredictDB training "~ R^2 ~ " in MESA_AFA")) +
        ylab(bquote("PrediXcan " ~ R^2 ~ " in SAGE"))

sage.mesa.afhi.r2 = r2.all %>%
    filter(Prediction.Weights %in% c("MESA_AFHI", "Test: MESA_AFHI")) %>%
    dcast(., Gene ~ Prediction.Weights) %>%
    filter(`Test: MESA_AFHI` > 0.2) %>%
    as.data.table %>%
    rename(., "Test_R2" = "Test: MESA_AFHI")

sage.mesa.afhi.plot = sage.mesa.afhi.r2 %>% 
    ggplot(., aes(x = Test_R2, y = MESA_AFHI)) +
        geom_point(color = "goldenrod") +
        # add regression line (blue) to compare points
        geom_smooth(se = TRUE, method = "lm") +
        # add y=x line (red dotted)
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        ggtitle("PredictDB performance in SAGE with MESA_AFHI weights") +
        xlab(bquote("PredictDB training "~ R^2 ~ " in MESA_AFHI")) +
        ylab(bquote("PrediXcan " ~ R^2 ~ " in SAGE"))

sage.mesa.cau.r2 = r2.all %>%
    filter(Prediction.Weights %in% c("MESA_CAU", "Test: MESA_CAU")) %>%
    dcast(., Gene ~ Prediction.Weights) %>%
    filter(`Test: MESA_CAU` > 0.2) %>%
    as.data.table %>%
    rename(., "Test_R2" = "Test: MESA_CAU")

sage.mesa.cau.plot = sage.mesa.cau.r2 %>% 
    ggplot(., aes(x = Test_R2, y = MESA_CAU)) +
        geom_point(color = "purple") +
        # add regression line (blue) to compare points
        geom_smooth(se = TRUE, method = "lm") +
        # add y=x line (red dotted)
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        ggtitle("PredictDB performance in SAGE with MESA_CAU weights") +
        xlab(bquote("PredictDB training "~ R^2 ~ " in MESA_CAU")) +
        ylab(bquote("PrediXcan " ~ R^2 ~ " in SAGE"))


# save plots to file
sage.gtex7.plot.path = file.path(output.dir, paste0("sage.predixcan.r2.gtex7.dotplot.", plot.type))
sage.mesa.all.plot.path = file.path(output.dir, paste0("sage.predixcan.r2.mesaall.dotplot.", plot.type))
sage.mesa.afa.plot.path = file.path(output.dir, paste0("sage.predixcan.r2.mesaafa.dotplot.", plot.type))
sage.mesa.afhi.plot.path = file.path(output.dir, paste0("sage.predixcan.r2.mesaafhi.dotplot.", plot.type))
sage.mesa.cau.plot.path = file.path(output.dir, paste0("sage.predixcan.r2.mesacau.dotplot.", plot.type))

ggsave(plot = sage.gtex7.plot, filename = sage.gtex7.plot.path, dpi = 300, type = "cairo") 
ggsave(plot = sage.mesa.all.plot, filename = sage.mesa.all.plot.path, dpi = 300, type = "cairo") 
ggsave(plot = sage.mesa.afa.plot, filename = sage.mesa.afa.plot.path, dpi = 300, type = "cairo") 
ggsave(plot = sage.mesa.afhi.plot, filename = sage.mesa.afhi.plot.path, dpi = 300, type = "cairo") 
ggsave(plot = sage.mesa.cau.plot, filename = sage.mesa.cau.plot.path, dpi = 300, type = "cairo") 

# save Rdata file
save.image(rdata.path)

# write results to output table
results.table.path = file.path(output.dir, "sage.predixcan.all.gene.results.txt")
fwrite(all.results, file = results.table.path, sep = "\t", quote = FALSE)
