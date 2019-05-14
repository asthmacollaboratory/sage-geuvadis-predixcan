library(ggplot2)
library(data.table)
library(dplyr)

# save a point jitter for (reproducibly) separating points on plots 
my.jitter = position_jitter(width = 0.1, height = 0.0001, seed = 2019)

# read results files
x = fread("twas.sims.all.results.2018-05-06.txt")

# should ideally confirm that these are true in results
nmodels = 11
nseeds  = 100
ngenes  = 98

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
my.colors = c(
    "AA to AA"   = "black",
    "CEU to AA"  = "blue",
    "YRI to AA"  = "red",
    "AA to CEU"  = "black",
    "CEU to CEU" = "blue",
    "YRI to CEU" = "red",
    "AA to YRI"  = "black",
    "CEU to YRI" = "blue",
    "YRI to YRI" = "red"
)

my.shapes = c(
    "AA to AA"   = "circle",
    "CEU to AA"  = "circle",
    "YRI to AA"  = "circle",
    "AA to CEU"  = "square",
    "CEU to CEU" = "square",
    "YRI to CEU" = "square",
    "AA to YRI"  = "triangle",
    "CEU to YRI" = "triangle",
    "YRI to YRI" = "triangle"
)

my.fillcolors = c("AA" = "gray95", "CEU" = "steelblue", "YRI" = "red")
my.linecolors = c("AA" = "black", "CEU" = "darkblue", "YRI" = "firebrick4")
my.linetypes.trainpop = c("AA" = "solid", "CEU" = "dashed", "YRI" = "dotted")

my.fillcolors.9 = c(
    "AA to AA"   = "gray95",
    "CEU to AA"  = "steelblue",
    "YRI to AA"  = "red",
    "AA to CEU"  = "gray95",
    "CEU to CEU" = "steelblue",
    "YRI to CEU" = "red",
    "AA to YRI"  = "gray95",
    "CEU to YRI" = "steelblue",
    "YRI to YRI" = "red"
)

my.linecolors.9 = c(
    "AA to AA"   = "black",
    "CEU to AA"  = "black",
    "YRI to AA"  = "black",
    "AA to CEU"  = "darkblue",
    "CEU to CEU" = "darkblue",
    "YRI to CEU" = "darkblue",
    "AA to YRI"  = "firebrick4",
    "CEU to YRI" = "firebrick4",
    "YRI to YRI" = "firebrick4"
)

my.breaks.diffpop = c(
    "YRI to AA",
    "AA to YRI",
    "AA to CEU",
    "CEU to AA",
    "YRI to CEU",
    "CEU to YRI"
)

my.breaks.allpop = c(
    "AA to AA",
    "CEU to CEU",
    "YRI to YRI",
    "YRI to AA",
    "AA to YRI",
    "AA to CEU",
    "CEU to AA",
    "YRI to CEU",
    "CEU to YRI"
)

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

# custom theme function
theme_klk = function(){
    theme_gray() +
    theme(
         text = element_text(size = 16),
         panel.background = element_rect(fill = "white"),
         panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "lightgrey"),
         panel.grid.minor = element_line(size = 0.25, linetype = "dashed", color = "lightgrey"),
         legend.position = "right"
    )
}


# boxplot of t-statistics for each train-test scenario
# serves as proxy for how p-values would behave
g1 = x %>%
    dplyr::filter(
        Original_Model == 0.1
    ) %>% 
    mutate(Train_Test = paste(Train_Pop, Test_Pop, sep = " to ")) %>% 
    mutate(
        Train_Test = factor(Train_Test, levels = my.breaks.allpop.reverse)
    ) %>%
    ggplot(., aes(x = Train_Test, y = T_value, fill = Train_Pop, color = Test_Pop))  +
        geom_boxplot(
            outlier.color = "black",
            lwd = 1 
        ) + 
        ggtitle(
            bquote("Distributions of "~t~"-statistics from TWAS association tests in AA, CEU, and YRI"),
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTLs per gene and effect size "~beta~" = 0.1")
        ) + 
        ylab(bquote(t~"-statistic")) +
        xlab("Train - Test scenario") +
        scale_fill_manual(values  = my.fillcolors) +
        scale_color_manual(values = my.linecolors) + 
       #  scale_linetype_manual(values = my.linetypes.trainpop) +
        theme_klk() +
        theme(
             legend.position = "none",
             axis.text.x = element_text(angle = 90, hjust = 1)
        ) +
        facet_grid(
            ~ Prop_Shared_eQTL,
            labeller = labeller(Prop_Shared_eQTL = c(
                `0` = "0% shared eQTL",
                `0.5` = "50% shared eQTL",
                `0.9` = "90% shared eQTL",
                `1`   = "100% shared eQTL"
                )
            )
        )
ggsave(g1, file = "twas.sim.results.2019-05-06.tstatistic.boxplot.png", units = "in", width = 18, height = 6)

g2 = x %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Power, P_value) %>%
    group_by(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL) %>%
    summarize(Post_Hoc_Power = sum(P_value < 0.05/98)) %>%
    dplyr::mutate(
        Train_Test = paste(Train_Pop, Test_Pop, sep = " to ")
    ) %>%
    dplyr::filter(
        Original_Model != 0
    ) %>%
    ggplot(., aes(x = log10(Original_Model), y = Post_Hoc_Power, color = Train_Test))  +
        geom_point(aes(shape = Train_Test), alpha = 0.2, position = my.jitter) +
        geom_smooth(
            aes(linetype = Train_Test),
            method = "glm",
            method.args = list(family = "binomial"),
            se = TRUE,
            alpha = 0.1,
            span = 1.0
        ) +
        scale_linetype_manual(
            name = "Train to Test",
            values = my.linetypes,
            breaks = my.breaks.allpop
        ) +
        scale_color_manual(
            name = "Train to Test",
            values = my.colors,
            breaks = my.breaks.allpop
        ) +
        scale_shape_manual(
            name = "Train to Test",
            values = my.shapes,
            breaks = my.breaks.allpop
        ) +
        facet_grid(~ Prop_Shared_eQTL,
            labeller = labeller(Prop_Shared_eQTL = c(
                `0` = "0% shared eQTL",
                `0.5` = "50% shared eQTL",
                `0.9` = "90% shared eQTL",
                `1`   = "100% shared eQTL"
                )
            )
        ) +
        theme_klk() +
        theme(legend.key.width = unit(1, "in")) + 
        guides(
            fill     = guide_legend(keywidth = 1, keyheight = 1),
            linetype = guide_legend(keywidth = 6, keyheight = 1),
            color    = guide_legend(keywidth = 6, keyheight = 1)
        ) +
        ylab("Power") +
        xlab(bquote(log[10]~"("~beta~")"))

g3 = x %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Power, P_value) %>%
    group_by(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL) %>%
    summarize(Post_Hoc_Power = sum(P_value < 0.05/98)) %>%
    dplyr::mutate(
        Train_Test = paste(Train_Pop, Test_Pop, sep = " to ")
    ) %>%
    dplyr::filter(
        Original_Model != 0,
        Train_Pop != Test_Pop,
        Prop_Shared_eQTL != 1
    ) %>%
    ggplot(., aes(x = log10(Original_Model), y = Post_Hoc_Power, color = Train_Test)) +
        geom_point(aes(shape = Train_Test), alpha = 0.2, position = my.jitter) +
        geom_smooth(aes(linetype = Train_Test), method = "glm", method.args = list(family = "binomial"), se = TRUE, alpha = 0.1, span = 0.5) +
        geom_vline(xintercept = log10(0.1), color = "red", linetype = "dotdash") +
        scale_shape_manual(
            name   = "Train to Test",
            values = my.shapes,
            breaks = my.breaks.diffpop
        ) +
        scale_linetype_manual(
            name   = "Train to Test",
            values = my.linetypes,
            breaks = my.breaks.diffpop
        )  +
        scale_color_manual(
            name   = "Train to Test",
            values = my.colors,
            breaks = my.breaks.diffpop
        ) +
        facet_grid(~ Prop_Shared_eQTL,
            labeller = labeller(Prop_Shared_eQTL = c(
                `0` = "0% shared eQTL",
                `0.5` = "50% shared eQTL",
                `0.9` = "90% shared eQTL",
                `1`   = "100% shared eQTL"
                )
            )
        ) +
        theme_klk() +
        theme(legend.key.width = unit(1, "in")) + 
        guides(
            fill     = guide_legend(keywidth = 1, keyheight = 1),
            linetype = guide_legend(keywidth = 6, keyheight = 1),
            color    = guide_legend(keywidth = 6, keyheight = 1)
        ) +
        ylab("Power") +
        xlab(bquote(log[10]~"("~beta~")")) +
        ggtitle(
            bquote("Power of TWAS association tests with cross-population predicted expression"),
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTLs per gene") 
        )
ggsave(g3, file = "twas.sim.results.2019-05-06.power.curves.png", units = "in", width = 18, height = 6)

g4 = x %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Power, P_value) %>%
    group_by(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL) %>%
    summarize(Post_Hoc_Power = sum(P_value < 0.05/98)) %>%
    dplyr::mutate(
        Train_Test = paste(Train_Pop, Test_Pop, sep = " to ")
    ) %>%
    dplyr::filter(
        Original_Model != 0,
    ) %>%
    ggplot(., aes(x = log10(Original_Model), y = Post_Hoc_Power, color = Train_Test)) +
        geom_point(aes(shape = Train_Test), alpha = 0.2, position = my.jitter) +
        geom_smooth(aes(linetype = Train_Test), method = "glm", method.args = list(family = "binomial"), se = TRUE, alpha = 0.1, span = 0.5) +
        geom_vline(xintercept = log10(0.1), color = "red", linetype = "dotdash") +
        scale_shape_manual(
            name   = "Train to Test",
            values = my.shapes,
            breaks = my.breaks.allpop
        ) +
        scale_linetype_manual(
            name   = "Train to Test",
            values = my.linetypes,
            breaks = my.breaks.allpop
        )  +
        scale_color_manual(
            name   = "Train to Test",
            values = my.colors,
            breaks = my.breaks.allpop
        ) +
        facet_grid(~ Prop_Shared_eQTL,
            labeller = labeller(Prop_Shared_eQTL = c(
                `0`   = "0% shared eQTL",
                `0.5` = "50% shared eQTL",
                `0.9` = "90% shared eQTL",
                `1`   = "100% shared eQTL"
                )
            )
        ) +
        theme_klk() +
        theme(legend.key.width = unit(1, "in")) + 
        guides(
            fill     = guide_legend(keywidth = 1, keyheight = 1),
            linetype = guide_legend(keywidth = 6, keyheight = 1),
            color    = guide_legend(keywidth = 6, keyheight = 1)
        ) +
        ylab("Power") +
        xlab(bquote(log[10]~"("~beta~")")) +
        ggtitle(
            bquote("Power of TWAS association tests with cross-population predicted expression"),
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTLs per gene") 
        )
ggsave(g4, file = "twas.sim.results.2019-05-06.power.curves.allpops.png", units = "in", width = 18, height = 6)


# at this juncture, it is handy to save another data.frame to summarize power calculations
# will use this for future plots
x.power = x %>%
    dplyr::filter(Original_Model == 0.1) %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Same_Causal_eQTL, P_value) %>%
    group_by(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Same_Causal_eQTL) %>%
    summarize(
        Post_Hoc_Power = sum(P_value < 0.05/ngenes) / (nseeds),
    ) %>%
    dplyr::mutate(
        Train_Test = paste(Train_Pop, Test_Pop, sep = " to ")
    ) %>%
    dplyr::mutate(
        Train_Test = factor(Train_Test, levels = my.breaks.allpop.reverse)
    ) %>%
    dplyr::filter(
        Original_Model != 0
    ) %>%
    ungroup %>%
    select(Train_Pop, Test_Pop, Train_Test, Prop_Shared_eQTL, Post_Hoc_Power) %>%
    group_by(Train_Test, Prop_Shared_eQTL) %>%
    summarize(
        Power = sum(Post_Hoc_Power),
        Power_SE = sd(Post_Hoc_Power) / sqrt(n())
    )

# busy bar chart, including all variables, useful for supplement
g5 = ggplot(x.power, aes(x = Train_Test, y = Power, fill = Train_Test, color = Train_Test, linetype = Train_Test))  +
        scale_fill_manual(values  = my.fillcolors.9) +
        scale_color_manual(values = my.linecolors.9) +
        scale_linetype_manual(values = my.linetypes) +
        geom_bar(stat = "identity", size = 2) +
        theme_klk() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        ) + 
        guides(
            fill     = guide_legend(keywidth = 1, keyheight = 1),
            linetype = guide_legend(keywidth = 6, keyheight = 1),
            color    = guide_legend(keywidth = 6, keyheight = 1)
        ) +
        ylab("Power") +
        xlab("Train-test scenario") +
        facet_grid(
            ~ Prop_Shared_eQTL,
            labeller = labeller(Prop_Shared_eQTL = c(
                `0` = "0% shared eQTL",
                `0.5` = "50% shared eQTL",
                `0.9` = "90% shared eQTL",
                `1`   = "100% shared eQTL"
                )
            )
        )

g6 = x.power %>%
    dplyr::filter(
        !(Train_Test %in% c("AA to AA", "CEU to CEU", "YRI to YRI")) &
        Prop_Shared_eQTL != 1
    ) %>%
    ungroup %>%
    mutate(
        Train_Test = factor(Train_Test, levels = c("CEU to YRI", "YRI to CEU", "CEU to AA", "AA to CEU", "AA to YRI", "YRI to AA")),
        Prop_Shared_eQTL = recode(Prop_Shared_eQTL, `0` = "0% shared eQTL", `0.5` = "50% shared eQTL", `0.9` = "90% shared eQTL")
    ) %>%
    ggplot(., aes(x = Train_Test, y = Power, group = Prop_Shared_eQTL, fill = Prop_Shared_eQTL))  +
        geom_bar(stat = "identity", position = position_dodge(), size = 2) +
        scale_color_manual(name = "Shared eQTL", values = c(`0% shared eQTL` = "lightblue", `50% shared eQTL` = "steelblue", `90% shared eQTL` = "darkblue"), aesthetics = c("fill")) +
        theme_klk() + 
        ylab("Power") +
        xlab("Train-test scenario") +
        ggtitle(
            bquote("Power of TWAS association tests in AA, CEU, and YRI"),
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTLs per gene and effect size "~beta~" = 0.1")
        )

ggsave(g6, file = "twas.sim.results.2019-05-06.power.barchart.png", units = "in", width = 18, height = 6)

# more complete copy of previous plot, perhaps useful for demonstration purposes
g7 = x.power %>%
    ungroup %>%
    mutate(
        Train_Test = factor(Train_Test, levels = c("CEU to YRI", "YRI to CEU", "CEU to AA", "AA to CEU", "AA to YRI", "YRI to AA", "AA to AA", "CEU to CEU", "YRI to YRI")),
        Prop_Shared_eQTL = recode(Prop_Shared_eQTL, `0` = "0% shared eQTL", `0.5` = "50% shared eQTL", `0.9` = "90% shared eQTL", `1` = "All shared eQTL")
    ) %>%
    ggplot(., aes(x = Train_Test, y = Power, group = Prop_Shared_eQTL, fill = Prop_Shared_eQTL))  +
        geom_bar(stat = "identity", position = position_dodge(), size = 2) +
        scale_color_manual(
            name = "Shared eQTL",
            values = c(`0% shared eQTL` = "lightblue", `50% shared eQTL` = "steelblue", `90% shared eQTL` = "darkblue", `All shared eQTL` = "black"),
            aesthetics = c("fill")
        ) +
        theme_klk() + 
        guides(
            fill     = guide_legend(keywidth = 1, keyheight = 1),
            linetype = guide_legend(keywidth = 6, keyheight = 1),
            color    = guide_legend(keywidth = 6, keyheight = 1)
        ) +
        ylab("Power") +
        xlab("Train-test scenario") +
        ggtitle(
            bquote("Power of TWAS association tests in AA, CEU, and YRI"),
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTLs per gene and effect size "~beta~" = 0.1")
        )

ggsave(g7, file = "twas.sim.results.2019-05-06.power.barchart.allpops.png", units = "in", width = 18, height = 6)

## OLD CODE ###


## make factor of train --> test scenarios
#twas.results = twas.results %>%
#    mutate(
#        Train_Test = factor(
#            paste(Train_Pop, Test_Pop, sep = " to "),
#            levels = c(
#                "YRI to AA",
#                "AA to YRI",
#                "AA to CEU",
#                "CEU to AA",
#                "YRI to CEU",
#                "CEU to YRI",
#                "AA to AA",
#                "CEU to CEU",
#                "YRI to YRI"
#            )
#        )
#    )
#
## facet plot of distributions of t-statistics
## 9 facets in 3x3 grid
## rows are Test_Pop, cols are Train_pop 
#g = twas.results %>%
#    ggplot(., aes(x = T_value)) +
#        geom_histogram(aes(y = ..density..), binwidth=0.2) +
#        geom_density() +
#        facet_grid(cols = vars(Train_Pop), rows = vars(Test_Pop)) +
#        xlim(-5, 5) +
#        ggtitle(bquote("Distribution of "~t~"-statistics from association tests in AA, CEU, and YRI"),
#            subtitle = "Expression variously imputed from AA, CEU, and YRI") +
#        geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
#        xlab(bquote(t~"-statistic")) +
#        ylab("Density")
#tplot.histogram.path = file.path(output.dir, paste(output.filepfx, "tstatistics", "histogram", plot.type, sep = ".")) 
#ggsave(plot = g, filename = tplot.histogram.path, type = "cairo")
##mutate(Train_Test = paste(Train_Pop, Test_Pop, sep = "_")) %>%
#
## boxplots of distributions of t-statistics
## rows are Test_Pop, cols are Train_pop 
#my.fillcolors = c("AA" = "gray95", "CEU" = "steelblue", "YRI" = "red")
#my.linecolors = c("AA" = "black", "CEU" = "darkblue", "YRI" = "firebrick4")
#my.linetypes = c("AA" = "solid", "CEU" = "dashed", "YRI" = "dotted")
#g = twas.results %>%
#    ggplot(., aes(
#            x = Train_Test,
#            y = T_value,
#            fill = Train_Pop,
#            color = Test_Pop)
#        ) +
#        geom_boxplot(
#            outlier.color = "black",
#            lwd = 1
#        ) +
#        ggtitle(
#            bquote("Distributions of "~t~"-statistics from association tests in AA, CEU, and YRI"),
#            subtitle = "Expression variously imputed from AA, CEU, and YRI"
#        ) +
#        ylab(bquote(t~"-statistic")) +
#        xlab("Train - Test scenario") +
#        scale_fill_manual(values  = my.fillcolors) +
#        scale_color_manual(values = my.linecolors) + 
#      #  scale_linetype_manual(values = my.linetypes) +
#        theme(
#            text = element_text(size = 16),
#            panel.background = element_rect(fill = "white"),
#            panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "lightgrey"),
#            panel.grid.minor = element_line(size = 0.25, linetype = "dashed", color = "lightgrey"),
#            legend.position = "none"
#        )
#
#
#tplot.boxplot.path = file.path(output.dir, paste(output.filepfx, "tstatistics", "boxplot", plot.type, sep = ".")) 
#ggsave(plot = g, filename = tplot.boxplot.path, type = "cairo", width = 10, height = 15, unit = "in")
#            #linetype = Test_Pop)
#            #color = Test_Pop)
#        #scale_fill_manual(values = my.colors) +
#        #geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
#
## boxplots of -log10 p-values
#g = twas.results %>%
#    ggplot(., aes(
#            x = Train_Test,
#            y = -log10(P_value),
#            fill = Train_Pop,
#            color = Test_Pop)
#        ) +
#        geom_boxplot(
#            outlier.color = "black",
#            lwd = 1 
#        ) +
#        ggtitle(
#            bquote("Distributions of "~p~"-values from association tests in AA, CEU, and YRI"),
#            subtitle = "Expression variously imputed from AA, CEU, and YRI"
#        ) +
#        ylab(bquote(-log10~"("~p~"-value)")) +
#        xlab("Train - Test scenario") +
#        scale_fill_manual(values  = my.fillcolors) +
#        scale_color_manual(values = my.linecolors) + 
#      #  scale_linetype_manual(values = my.linetypes) +
#        theme(
#            text = element_text(size = 16),
#            panel.background = element_rect(fill = "white"),
#            panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "lightgrey"),
#            panel.grid.minor = element_line(size = 0.25, linetype = "dashed", color = "lightgrey"),
#            legend.position = "none"
#        ) +
#        geom_hline(yintercept = -log10(0.05/98), color = "red", linetype = "twodash") +
#        geom_point(
#            data = subset(twas.results, Original_Model != 0),
#            aes(x = Train_Test, y = -log10(P_value)),
#            color = "purple", size = 5, shape = 15
#        ) 
#tplot.pvalue.path = file.path(output.dir, paste(output.filepfx, "pvalues", "boxplot", plot.type, sep = ".")) 
#ggsave(plot = g, filename = tplot.pvalue.path, type = "cairo", width = 10, height = 15, unit = "in")
#
## save a version that better shows the p-value distribution without outliers
#g = g + ylim(0,6)
#tplot.pvalue.path = file.path(output.dir, paste(output.filepfx, "pvalues", "boxplot", "zoomed", plot.type, sep = ".")) 
#ggsave(plot = g, filename = tplot.pvalue.path, type = "cairo", width = 10, height = 15, unit = "in")
#
#
## QQ plot of p-values, also in facet
#
## make a transformed copy of the p-values for plotting
#xlim = NULL
#ylim = NULL
#df.copy = data.frame(
#    #"expected"  = -log10(ppoints(length(twas.results$P_value))),
#    #"expected"  = -log10(ppoints(rep(c(1:ngenes), 9))),
#    "P_value"  = twas.results$P_value,
#    "Train_Pop" = twas.results$Train_Pop,
#    "Test_Pop"  = twas.results$Test_Pop
#)
#df.copy = df.copy %>%
#    group_by(Train_Pop, Test_Pop) %>%
#    arrange(Train_Pop, Test_Pop, desc(P_value)) %>%
#    mutate(
#        "observed" = -log10(P_value),
#        "expected" = -log10(ppoints(ngenes))
#    )
#
#g = ggplot(df.copy, aes(x = expected, y = observed)) + geom_point()
#
## draw the qq-line
#g = g + geom_abline(intercept = 0, slope = 1, color = "red")
#
## prettify the axis labels
#g = g + xlab("Theoretical Quantiles") + ylab("Sample Quantiles")
#
## adjust the axis limits
#g = g + coord_cartesian(xlim = xlim, ylim = ylim)
#
## apply facet
#g = g + facet_grid(cols = vars(Train_Pop), rows = vars(Test_Pop))
#
#qqplot.path = file.path(output.dir, paste(output.filepfx, "qqplot", plot.type, sep = ".")) 
#ggsave(plot = g, filename = qqplot.path, type = "cairo") 

