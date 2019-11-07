##!/usr/bin/env Rscript --vanilla
# ==========================================================================================
# coded by Kevin L. Keys (2018)
# ==========================================================================================


# ==========================================================================================
# libraries
# ==========================================================================================
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(methods))
suppressMessages(library(optparse))
suppressMessages(library(dunn.test))

# parse command line variables
option_list = list(
    make_option(
        c("-r", "--results-file"),
        type    = "character",
        default = NULL,
        help    = "The data frame containing all results to analyze.",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-directory"),
        type    = "character",
        default = NULL,
        help    = "Prefix for output files.",
        metavar = "character"
    ),
    make_option(
        c("-p", "--plot-filetype"),
        type    = "character",
        default = "pdf",
        help    = "Output filetype for the plots, e.g. 'pdf', 'png', etc. [default: %default]",
        metavar = "character"
    ),
    make_option(
        c("-k", "--num-eQTL"),
        type    = "integer",
        default = "10",
        help    = "Number of causal eQTL used in the prediction models [default: %default]",
        metavar = "character"
    )
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

# ==========================================================================================
# script variables
# ==========================================================================================
results.file  = opt$results_file
output.dir    = opt$output_directory
plot.filetype = opt$plot_filetype
k.value = opt$num_eQTL


# plots for k == 1 must be produced separately
# use variable K to control which plot is produced
K = c(10,20,40)

# ==========================================================================================
# subroutines
# ==========================================================================================

# subroutine to make faceted corr plot
make.corr.plot = function(x, k.val = 1, binwidth = 0.01, xlim.lo = -1, xlim.hi = 1, ylim.lo = 0, ylim.hi = 25000){
    corr.plot = x %>%
        dplyr::filter(k == k.val) %>%
        ggplot(aes(x = Correlation)) +
            geom_histogram(binwidth = binwidth, color = "black", fill = "white") +
            ggtitle("Cross-population imputation accuracy", subtitle = paste0("Correlation for causal number of eQTLs k = ", k.val)) +
            xlab(expression(Spearman~italic(rho))) +
            xlim(xlim.lo, xlim.hi) +
            ylim(ylim.lo, ylim.hi)
    return(corr.plot)
}

# subroutine to make faceted R2 plot
make.r2.plot = function(x, k.val = 1, binwidth = 0.005, xlim.lo = 0, xlim.hi = 1, ylim.lo = 0, ylim.hi = 25000){
    r2.plot = x %>%
        dplyr::filter(k == k.val) %>%
        ggplot(aes(x = R2)) +
            geom_histogram(binwidth = binwidth, color = "black", fill = "white") +
            ggtitle("Cross-population imputation accuracy", subtitle = paste0("R2 for causal number of eQTLs k = ", k.val)) +
            xlab(expression(italic(R)^2)) +
            xlim(xlim.lo, xlim.hi) +
            ylim(ylim.lo, ylim.hi)
    return(r2.plot)
}



# ==========================================================================================
# script code
# ==========================================================================================

# fix various plotting parameters here
# writing these once and reusing makes for tidier code
my.linetypes = c(
    "AA to AA"   = "solid",
    "CEU to AA"  = "solid",
    "YRI to AA"  = "solid",
    "AA to CEU"  = "dashed",
    "CEU to CEU" = "dashed",
    "YRI to CEU" = "dashed",
    "AA to YRI"  = "dotted",
    "CEU to YRI" = "dotted",
    "YRI to YRI" = "dotted"
)

my.breaks.allpop = c(
    "YRI to AA",
    "AA to YRI",
    "AA to CEU",
    "CEU to AA",
    "YRI to CEU",
    "CEU to YRI",
    "AA to AA",
    "CEU to CEU",
    "YRI to YRI"
)

my.colors.withgrey = c(
    "AA to AA"   = "darkgrey",
    "CEU to AA"  = "blue",
    "YRI to AA"  = "red",
    "AA to CEU"  = "black",
    "CEU to CEU" = "darkgrey",
    "YRI to CEU" = "red",
    "AA to YRI"  = "black",
    "CEU to YRI" = "blue",
    "YRI to YRI" = "darkgrey"
)

# enable plotting of PNG and similar files
options(bitmapType = 'cairo')

# load data
x = fread(results.file)

# add a train-test column
x$Train_Test = paste(x$Train_Pop, x$Test_Pop, sep = " to ")

# make plots for each k
# not sure of most efficient way to do this
# maybe make list of ggplots?
# copy boilerplate code for now
if (k.value %in% K) {

    # filepaths for saving plots
    r2.facet.plot.path          = file.path(output.dir, paste0("1kg.sims.corrs.facet.k", k.value, ".", plot.filetype))
    corr.facet.plot.path        = file.path(output.dir, paste0("1kg.sims.R2.facet.k", k.value, ".", plot.filetype))
    corr.by.propsharedeqtl.path = file.path(output.dir, paste0("1kg.sims.corrs.by.propsharedeqtl.k", k.value, ".", plot.filetype))

    # make facet plots
    r2.facet.plot   = make.r2.plot(x, k.val = k.value)
    corr.facet.plot = make.corr.plot(x, k.val = k.value)

    # can only plot regression lines for one value of k at a time
    corr.by.propsharedeqtl = x %>%
        dplyr::filter(
            (k == k.value) &
            (CEU_prop == 0.2) &
            (YRI_prop == 0.8)
            #(Train_Pop != Test_Pop) &
        ) %>%
        group_by(Train_Test, prop_shared_eqtl) %>%
        summarize(Corr_Mean = mean(Correlation, na.rm = TRUE), Corr_StdErr = sd(Correlation, na.rm = TRUE)) %>%
        as.data.table

    # make regression lines of correlation by prop_shared_eqtl
    corr.by.propsharedeqtl.plot = ggplot(corr.by.propsharedeqtl, aes(x = prop_shared_eqtl, y = Corr_Mean, group = Train_Test, color = Train_Test)) +
        geom_point(alpha = 0.05) +
        geom_smooth(aes(linetype = Train_Test), se = TRUE, method = "loess", size = 1.0, level = 0.95) +
        xlab("Proportion of shared eQTLs") +
        ylab("Spearman Correlation") +
        ggtitle("Crosspopulation correlations of predicted versus simulated gene expression", subtitle = paste0("Number of causal eQTLs: ", k.value)) +
        xlim(0, 1) +
        ylim(0, 0.4) +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "lightgrey"),
            panel.grid.minor = element_line(size = 0.25, linetype = "dashed", color = "lightgrey")
        ) +
        scale_linetype_manual(
            name   = "Train to Test",
            values = my.linetypes, 
            breaks = my.breaks.allpop 
        )  +
        scale_color_manual(
            name   = "Train to Test",
            values = my.colors.withgrey, 
            breaks = my.breaks.allpop  
        ) +
        guides(
            fill     = guide_legend(keywidth = 1, keyheight = 1),
            linetype = guide_legend(keywidth = 6, keyheight = 1),
            color    = guide_legend(keywidth = 6, keyheight = 1)
        ) +
        geom_hline(yintercept = sqrt(0.15), size = 1.5, color = "black", linetype = "solid")

    # save all plots to file
    ggsave(filename = corr.facet.plot.path, plot = corr.facet.plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(filename = r2.facet.plot.path, plot = r2.facet.plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(filename = corr.by.propsharedeqtl.path, plot = corr.by.propsharedeqtl.plot, dpi = 300, width = 20, height = 10, units = "in")
} else {
    # set factor levels for Train-Test labels
    # will use this to sort left-right order of bar chart
    my.breaks.allpop.reverse = c(
        "CEU to YRI",
        "YRI to CEU",
        "AA to CEU",
        "CEU to AA",
        "AA to YRI",
        "YRI to AA",
        "AA to AA",
        "CEU to CEU",
        "YRI to YRI"
    )

    # set file path to save plots
    corr.by.propsharedeqtl.path = file.path(output.dir, paste0("1kg.sims.corrs.by.propsharedeqtl.k", k.value, ".", plot.filetype))

    # produce *bar chart* for k == 1
    # this will cover the TRUE/FALSE nature of having 1 shared eQTL
    corr.by.propsharedeqtl = x %>%
        dplyr::filter(k == 1 & prop_shared_eqtl %in% c(0,1)) %>%
        group_by(Train_Pop, Test_Pop, gene, same_eqtls) %>%
        summarize(m = mean(Correlation, na.rm=T)) %>%
        ungroup %>%
        group_by(Train_Pop, Test_Pop, same_eqtls) %>%
        summarize(Correlation_Mean = mean(m), Correlation_StdErr = sd(m)) %>%
        as.data.table

    corr.by.propsharedeqtl.plot = corr.by.propsharedeqtl %>%
        mutate(
            Train_Test = factor(
                paste(Train_Pop, Test_Pop, sep = " to "),
                levels = my.breaks.allpop.reverse
            )
        ) %>%
        rename("Same eQTL" = same_eqtls) %>%
        ggplot(., aes(x = Train_Test, y = as.numeric(Correlation_Mean), fill = `Same eQTL`)) +
            geom_bar(position = position_dodge(), stat = "identity") +
            geom_errorbar(
                aes(ymin = Correlation_Mean - 1.96*Correlation_StdErr, ymax = Correlation_Mean + 1.96*Correlation_StdErr),
                position = position_dodge(0.9),
                size = 0.3,
                width = 0.2
            ) +
            ggtitle(
                "Crosspopulation correlations of predicted versus simulated gene expression",
                subtitle = paste0("Number of causal eQTL: ", k.value)) +
            theme(
                text = element_text(size = 30),
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "lightgrey"),
                panel.grid.minor = element_line(size = 0.25, linetype = "dashed", color = "lightgrey"),
                axis.text.x = element_text(angle = 90, hjust = 1)
            ) +
            xlab("Train-test scenario") +
            ylab("Spearman Correlation") +
            scale_fill_manual(
                name   = "Same eQTL",
                values = c("steelblue", "firebrick"),
                breaks = c("FALSE", "TRUE")
            )

    # save plot to file
    ggsave(filename = corr.by.propsharedeqtl.path, plot = corr.by.propsharedeqtl.plot, dpi = 300, width = 20, height = 10, units = "in")
}


# add special analyses and plots for k.value = 10
# the analyses will output to the terminal
if (k.value == 10) {

    # need some comparison statistics
    # are groups significantly different from each other?
    x.kruskal = x %>%
        # pull results for current k, purge pop-to-itself, and 0.90 < prop_shared_eqtl < 1.0
        dplyr::filter(
            (k == k.value) &
            (Train_Pop != Test_Pop) &
            (CEU_prop == 0.2) &
            (YRI_prop == 0.8)
        ) %>%
        kruskal.test(Correlation ~ as.factor(Train_Test), data = .)
    print(x.kruskal)

    # *which* groups are different from each other?
    x.dunn = x %>%
        # pull results for current k, purge pop-to-itself, and 0.90 < prop_shared_eqtl < 1.0
        dplyr::filter(
            (k == k.value) &
            (Train_Pop != Test_Pop) &
            (CEU_prop == 0.2) &
            (YRI_prop == 0.8)
        ) %>%
        group_by(Train_Test, gene, prop_shared_eqtl) %>%
        summarize(Correlation = mean(Correlation, na.rm = TRUE)) %>%
        select(Train_Test, gene, prop_shared_eqtl, Correlation) %>%
        dunn.test(x = .$Correlation, g = .$Train_Test, method = "bonferroni")
    dunn.results = data.table(data.frame(x.dunn[-1]))
    dunn.results.path = file.path(output.dir, paste0("1kg.sims.corrs.dunn.k", k.value, ".txt"))
    fwrite(x = dunn.results, file = dunn.results.path, sep = "\t")

    # want to produce table of results for imputing into same pop
    # this is useful as supp table for manuscript
    samepop.results = x %>%
        dplyr::filter(
            (Train_Pop == Test_Pop) &
            (CEU_prop == 0.2) &
            (YRI_prop == 0.8)
        ) %>%
        group_by(Train_Pop, Test_Pop, k, gene, prop_shared_eqtl) %>%
        summarize(Correlation = mean(Correlation, na.rm = TRUE)) %>%
        group_by(Train_Pop, Test_Pop, k, prop_shared_eqtl) %>%
        summarize(Correlation_Mean = mean(Correlation), Correlation_StdErr = sd(Correlation)) %>%
        as.data.table
    samepop.results.path = file.path(output.dir, paste0("1kg.sims.corrs.by.propsharedeqtl.samepop.summary.txt"))
    fwrite(x = samepop.results, file = samepop.results.path, sep = "\t")
}
