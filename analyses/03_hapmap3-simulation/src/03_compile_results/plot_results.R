##!/usr/bin/env Rscript --vanilla
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
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
        help    = "Output filetype for the plots, e.g. 'pdf', 'png', 'svg', etc.",
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


# which model sizes did we test?
# note that this is *FIXED* here; the script *WILL BREAK* if $K doesn't match prediction output!
K = c(1,5,10,20)

# ==========================================================================================
# subroutines
# ==========================================================================================

# subroutine to make faceted corr plot
make.corr.plot = function(x, k.val = 1, binwidth = 0.01, xlim.lo = -1, xlim.hi = 1, ylim.lo = 0, ylim.hi = 25000){
corr.plot = x %>%
    dplyr::filter(k == k.val) %>%
    ggplot(aes(x = Correlation)) +
        geom_histogram(binwidth = binwidth, color = "black", fill = "white") +
        facet_grid(Train_Pop ~ Test_Pop) + 
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
            facet_grid(Train_Pop ~ Test_Pop) +
            ggtitle("Cross-population imputation accuracy", subtitle = paste0("R2 for causal number of eQTLs k = ", k.val)) +
            xlab(expression(italic(R)^2)) +
            xlim(xlim.lo, xlim.hi) +
            ylim(ylim.lo, ylim.hi)
    return(r2.plot)
}



# ==========================================================================================
# script code 
# ==========================================================================================

# load data
x = fread(results.file)

# add a train-test column
x$Train_Test = paste(x$Train_Pop, x$Test_Pop, sep = " to ")

# make plots for each k
# not sure of most efficient way to do this
# maybe make list of ggplots?
# copy boilerplate code for now
for (k.val in K) {

    # filepaths for saving plots
    r2.facet.plot.path          = file.path(output.dir, paste0("1kg.sims.corrs.facet.k", k.val, ".", plot.filetype)) 
    corr.facet.plot.path        = file.path(output.dir, paste0("1kg.sims.R2.facet.k", k.val, ".", plot.filetype))
    corr.by.propsharedeqtl.path = file.path(output.dir, paste0("1kg.sims.corrs.by.propsharedeqtl.k", k.val, ".", plot.filetype))

    # make facet plots
    r2.facet.plot   = make.r2.plot(x, k.val = k.val)
    corr.facet.plot = make.corr.plot(x, k.val = k.val)

    # make regression lines of correlation by prop_shared_eqtl
    corr.by.propsharedeqtl.plot = x %>%
        dplyr::filter(
            (Train_Pop != Test_Pop) &
            (k == k.val) &
            (prop_shared_eqtl < 0.91) 
        ) %>%
        ggplot(., aes(x=prop_shared_eqtl, y = Correlation, group = Train_Test, color = Train_Test)) +
            geom_smooth(aes(linetype = Train_Test), se = T, method = "lm", size = 2.5) +
            xlab("Proportion of shared eQTLs") +
            ylab("Spearman Correlation") +
            ggtitle("Crosspopulation correlations of predicted versus simulated gene expression", subtitle = paste0("Number of causal eQTLs: ", k.val)) +
            xlim(0, 0.9) +
            theme(
                text = element_text(size = 30),
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "lightgrey"),
                panel.grid.minor = element_line(size = 0.25, linetype = "dashed", color = "lightgrey")
            ) +
            scale_linetype_manual(
                name   = "Train to Test",
                values = c("dashed", "dotted", "solid", "dotted", "solid", "dashed"),
                breaks = c("YRI to AA", "AA to YRI", "AA to CEU", "CEU to AA", "YRI to CEU", "CEU to YRI")
            )  +
            scale_color_manual(
                name   = "Train to Test",
                values = c("black", "black", "blue", "blue", "red", "red"),
                breaks = c("YRI to AA", "AA to YRI", "AA to CEU", "CEU to AA", "YRI to CEU", "CEU to YRI")
            ) +
            guides(
                fill     = guide_legend(keywidth = 1, keyheight = 1),
                linetype = guide_legend(keywidth = 6, keyheight = 1),
                color    = guide_legend(keywidth = 6, keyheight = 1)
            )

    # save all plots to file
    ggsave(filename = corr.facet.plot.path, plot = corr.facet.plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(filename = r2.facet.plot.path, plot = r2.facet.plot, dpi = 300, width = 10, height = 10, units = "in")
    ggsave(filename = corr.by.propsharedeqtl.path, plot = corr.by.propsharedeqtl.plot, dpi = 300, width = 20, height = 10, units = "in")
}



# analyze just k = 10 for now
# the analyses will output to the terminal
k.val = 10

# need some comparison statistics
# are groups significantly different from each other? 
x.kruskal = x %>% 
	# pull results for current k, purge pop-to-itself, and 0.90 < prop_shared_eqtl < 1.0 
	dplyr::filter(
		(k == k.val) &
		(prop_shared_eqtl < 0.91) &
		(Train_Pop != Test_Pop)
	) %>%
	kruskal.test(Correlation ~ as.factor(Train_Test), data = .)
print(x.kruskal)

# *which* groups are different from each other?
x.dunn = x %>% 
	# pull results for current k, purge pop-to-itself, and 0.90 < prop_shared_eqtl < 1.0 
	dplyr::filter(
		(k == k.val) &
		(prop_shared_eqtl < 0.91) &  
		(Train_Pop != Test_Pop)
	) %>%
	dunn.test(x = .$Correlation, g = .$Train_Test, method = "bonferroni")
print(data.frame(x.dunn[-1]))
