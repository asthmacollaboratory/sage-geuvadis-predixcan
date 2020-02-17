suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(methods))
suppressMessages(library(optparse))

# parse command line variables
option_list = list(
    make_option(
        c("-a", "--results-file"),
        type    = "character",
        default = NULL, 
        help    = "The path to the results file, stored in tabular format", 
        metavar = "character"
    ),
    make_option(
        c("-b", "--output-directory"),
        type    = "character",
        default = NULL, 
        help    = "Directory where output files will be stored.", 
        metavar = "character"
    ),
    make_option(
        c("-c", "--plot-filetype"),
        type    = "character",
        default = "png", 
        help    = "Filetype for producing plots [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-d", "--random-seed"),
        type    = "integer",
        default = 2019, 
        help    = "Random seed to create reproducible plot jitters [default = %default].",
        metavar = "integer"
    ),
    make_option(
        c("-e", "--num-genes"),
        type    = "integer",
        default = 98, 
        help    = "Number of genes in the TWAS association tests [default = %default].",
        metavar = "integer"
    ),
    make_option(
        c("-f", "--num-models"),
        type    = "integer",
        default = 12, 
        help    = "Number of different causal effect sizes tested in TWAS association tests [default = %default].",
        metavar = "integer"
    ),
    make_option(
        c("-g", "--num-seeds"),
        type    = "integer",
        default = 100, 
        help    = "Number of different random instantiations in the TWAS association tests [default = %default].",
        metavar = "integer"
    ),
    make_option(
        c("-i", "--num-eQTL"),
        type    = "integer",
        default = 10, 
        help    = "Number of eQTL used to predict gene expression levels for TWAS association tests [default = %default].",
        metavar = "integer"
    ), ## note: this argument merely allows us to shoehorn this script through a unified QSUB submission interface
    make_option(
        c("-j", "--effect-size-to-plot"),
        type    = "numeric",
        default = 1e-2, 
        help    = "Effect size to use when making cross-sections of heritability[default = %default].",
        metavar = "numeric"
    ) ## note: this argument merely allows us to shoehorn this script through a unified QSUB submission interface
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

results.file = opt$results_file
output.dir   = opt$output_directory
plot.type    = opt$plot_filetype
seed         = opt$random_seed
effect.size  = opt$effect_size_to_plot

# should ideally confirm that these are true in results
nmodels      = opt$num_models
ngenes       = opt$num_genes 
nseeds       = opt$num_seeds 

# save a point jitter for (reproducibly) separating points on plots 
my.jitter = position_jitter(width = 0.1, height = 0.0001, seed = seed)

# read results files
x = fread(results.file)

# compute average genetic heritability for desired effect size to show
h2.effect.size = x %>%
    mutate(Train_Test = paste(Train_Pop, Test_Pop, sep = " to ")) %>% 
    group_by(Train_Test, Prop_Shared_eQTL, Original_Model) %>%
    summarize(h = mean(Heritability, na.rm = TRUE)) %>%
    dplyr::filter(Original_Model == effect.size) %>%
    ungroup %>%
    dplyr::select(h) %>%
    summarize(h2 = mean(h)) %>%
    as.numeric %>%
    round(., 3)

admix.effect.sizes = c(0.005, 0.01, 0.025)

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
cat("plotting g1\n\n")
g1 = x %>%
    dplyr::filter(
        Original_Model == effect.size &
        YRI_proportion == 0.8
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
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTL per gene and heritability "~h^2~" = "~.(h2.effect.size))
        ) + 
        ylab(bquote(t~"-statistic")) +
        xlab("Train - Test scenario") +
        scale_fill_manual(values  = my.fillcolors) +
        scale_color_manual(values = my.linecolors) + 
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
g1.path = file.path(output.dir, paste("twas.sim.results.tstatistic.boxplot", plot.type, sep = "."))
ggsave(g1, file = g1.path, units = "in", width = 18, height = 6)

g2 = x %>%
    dplyr::filter(
        Original_Model != 0 &
        YRI_proportion == 0.8
    ) %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Power, P_value) %>%
    group_by(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL) %>%
    summarize(Post_Hoc_Power = sum(P_value < 0.05/98)) %>%
    dplyr::mutate(
        Train_Test = paste(Train_Pop, Test_Pop, sep = " to ")
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
    dplyr::filter(
        Original_Model != 0 &
        Train_Pop != Test_Pop &
        Prop_Shared_eQTL != 1 &
        YRI_proportion == 0.8
    ) %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Power, P_value) %>%
    group_by(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL) %>%
    summarize(Post_Hoc_Power = sum(P_value < 0.05/98)) %>%
    dplyr::mutate(
        Train_Test = paste(Train_Pop, Test_Pop, sep = " to ")
    ) %>%
    ggplot(., aes(x = log10(Original_Model), y = Post_Hoc_Power, color = Train_Test)) +
        geom_point(aes(shape = Train_Test), alpha = 0.2, position = my.jitter) +
        geom_smooth(aes(linetype = Train_Test), method = "glm", method.args = list(family = "binomial"), se = TRUE, alpha = 0.1, span = 0.5) +
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
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTL per gene") 
        )
g3.path = file.path(output.dir, paste("twas.sim.results.power.curves", plot.type, sep = "."))
ggsave(g3, file = g3.path, units = "in", width = 18, height = 6)

# this applies a light grey color to pop-into-itself scenarios
my.colors.withgrey = c(
    "AA to AA"   = "lightgrey",
    "CEU to AA"  = "blue",
    "YRI to AA"  = "red",
    "AA to CEU"  = "black",
    "CEU to CEU" = "lightgrey",
    "YRI to CEU" = "red",
    "AA to YRI"  = "black",
    "CEU to YRI" = "blue",
    "YRI to YRI" = "lightgrey"
)

cat("plotting g4\n\n")
g4 = x %>%
    dplyr::filter(
        Original_Model != 0 &
        YRI_proportion == 0.8
    ) %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Power, P_value, Heritability) %>%
    group_by(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL) %>%
    summarize(
        Post_Hoc_Power = sum(P_value < 0.05/ngenes),
        h2 = mean(Heritability, na.rm = TRUE)
    ) %>%
    dplyr::mutate(
        Train_Test = paste(Train_Pop, Test_Pop, sep = " to ")
    ) %>%
    ggplot(., aes(x = h2, y = Post_Hoc_Power, color = Train_Test)) +
        geom_smooth(aes(linetype = Train_Test), method = "glm", method.args = list(family = "binomial"), se = TRUE, alpha = 0.1, span = 0.5) +
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
            values = my.colors.withgrey,
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
        xlab(bquote("Proportion of variance explained")) +
        #xlab(bquote("Heritability ("~h^2~")")) +
        ggtitle(
            bquote("Power of TWAS association tests with cross-population predicted expression"),
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTL per gene") 
        )
g4.path = file.path(output.dir, paste("twas.sim.results.power.curves.allpops", plot.type, sep = "."))
ggsave(g4, file = g4.path, units = "in", width = 18, height = 6)


# at this juncture, it is handy to save another data.frame to summarize power calculations
# will use this for future plots
x.power = x %>%
    dplyr::filter(
        Original_Model == effect.size &
        YRI_proportion == 0.8
    ) %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Same_Causal_eQTL, P_value, Heritability) %>%
    group_by(Train_Pop, Test_Pop, Seed, Original_Model, Prop_Shared_eQTL, Same_Causal_eQTL) %>%
    summarize(
        Post_Hoc_Power = sum(P_value < 0.05/ngenes) / (nseeds),
        h2 = mean(Heritability, na.rm = TRUE)
    ) %>%
    dplyr::mutate(
        Train_Test = factor(paste(Train_Pop, Test_Pop, sep = " to "), levels = my.breaks.allpop.reverse)
    ) %>%
    ungroup %>%
    select(Train_Test, Prop_Shared_eQTL, Post_Hoc_Power, h2) %>%
    group_by(Train_Test, Prop_Shared_eQTL) %>%
    summarize(
        Power = sum(Post_Hoc_Power),
        Power_SE = sd(Post_Hoc_Power) / sqrt(n()),
        Heritability_Mean = mean(h2, na.rm = TRUE)
    )

# busy bar chart, including all variables, useful for supplement
cat("plotting g5\n\n")
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

cat("plotting g6\n\n")
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
        geom_errorbar(
            aes(ymin = Power - 1.96*Power_SE, ymax = Power + 1.96*Power_SE),
            position = position_dodge(0.9),
            size = 0.3,
            width = 0.2 
        ) + 
        scale_color_manual(name = "Shared eQTL", values = c(`0% shared eQTL` = "lightblue", `50% shared eQTL` = "steelblue", `100% shared eQTL` = "darkblue"), aesthetics = c("fill")) +
        theme_klk() + 
        ylab("Power") +
        xlab("Train-test scenario") +
        ggtitle(
            bquote("Power of TWAS association tests in AA, CEU, and YRI"),
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTL per gene and "~h^2~" = "~.(h2.effect.size))
        )

g6.path = file.path(output.dir, paste("twas.sim.results.power.barcharts", plot.type, sep = "."))
ggsave(g6, file = g6.path, units = "in", width = 18, height = 6)

# more complete copy of previous plot, perhaps useful for demonstration purposes
cat("plotting g7\n\n")
g7 = x.power %>%
    ungroup %>%
    mutate(
        Train_Test = factor(Train_Test, levels = c("CEU to YRI", "YRI to CEU", "AA to CEU", "AA to YRI", "CEU to AA", "YRI to AA", "AA to AA", "CEU to CEU", "YRI to YRI")),
        Prop_Shared_eQTL = recode(Prop_Shared_eQTL, `0` = "0% shared eQTL", `0.5` = "50% shared eQTL", `0.9` = "90% shared eQTL", `1` = "All shared eQTL")
        #Train_Test = factor(Train_Test, levels = c("CEU to YRI", "YRI to CEU", "CEU to AA", "AA to CEU", "AA to YRI", "YRI to AA", "AA to AA", "CEU to CEU", "YRI to YRI")),
    ) %>%
    ggplot(., aes(x = Train_Test, y = Power, group = Prop_Shared_eQTL, fill = Prop_Shared_eQTL))  +
        geom_bar(stat = "identity", position = position_dodge(), size = 2) +
        scale_color_manual(
            name = "Shared eQTL",
            values = c(`0% shared eQTL` = "lightblue", `50% shared eQTL` = "steelblue", `All shared eQTL` = "darkblue"),
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
            subtitle = bquote("Expression imputed from AA, CEU, and YRI for "~k~" = 10 eQTL per gene and heritability "~h^2~" = "~.(h2.effect.size))
        )

g7.path = file.path(output.dir, paste("twas.sim.results.power.barcharts.allpops", plot.type, sep = "."))
ggsave(g7, file = g7.path, units = "in", width = 18, height = 6)

# produce power plots for varying admixture levels
x.power.admix = x %>%
    dplyr::filter(
        Original_Model %in% admix.effect.sizes &
        Prop_Shared_eQTL == 0.5 &
        Train_Pop != Test_Pop
    ) %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, YRI_proportion, Heritability, P_value) %>%
    mutate(
        Train_Test = factor(paste(Train_Pop, Test_Pop, sep = " to "), levels = my.breaks.allpop.reverse)
    ) %>%
    group_by(Train_Test, Seed, Original_Model, YRI_proportion) %>%
    summarize(
        Post_Hoc_Power = sum(P_value < 0.05/ngenes),
        h2 = mean(Heritability)
    ) %>%
    ungroup %>%
    group_by(Train_Test, Original_Model, YRI_proportion) %>%
    summarize(
        Power = sum(Post_Hoc_Power) / nseeds,
        Power_SE = sd(Post_Hoc_Power) / sqrt(n()),
        h2_Mean = mean(h2),
        h2_SE = sd(h2) / sqrt(n())
    ) %>%
    as.data.table

x.power.admix = x.power.admix %>%
    group_by(Original_Model) %>%
    summarize(Heritability_Mean = paste0("Heritability = ", round(mean(h2_Mean), 2))) %>%
    merge(., x.power.admix, by = "Original_Model")

cat("plotting g8\n\n")
g8 = ggplot(x.power.admix, aes(x = Train_Test, y = YRI_proportion, fill = Power)) +
    geom_tile() +
    geom_text(aes(label = round(Power, 2)), color = "white") +
    facet_grid(~ Heritability_Mean) +
    ggtitle("Power to detect causal genes in TWAS varies by admixture proportion",
        subtitle = bquote("Results for"~k~"= 10 eQTL, 50% eQTL shared across populations")
    ) +
    xlab("Train-test scenario") +
    ylab("YRI proportion in AA")

g8.path = file.path(output.dir, paste("twas.sim.results.power.heatmap.admixvary", plot.type, sep = "."))
ggsave(g8, file = g8.path, units = "in", width = 18, height = 6)


# as supplement to previous heatmap (g8), make a line plot with similar info
# this plot has added benefit of deemphasizing CEU <--> YRI relationships, which don't vary with admixture proportion 
x.power2 = x %>%
    dplyr::filter(
        Original_Model %in% admix.effect.sizes &
        Prop_Shared_eQTL == 0.5 &
        Train_Pop != Test_Pop
    ) %>%
    select(Train_Pop, Test_Pop, Seed, Original_Model, YRI_proportion, Heritability, P_value) %>%
    mutate(
        Train_Test = factor(paste(Train_Pop, Test_Pop, sep = " to "), levels = my.breaks.allpop.reverse)
    ) %>%
    mutate(
        Train_Test = dplyr::recode(Train_Test,
            `YRI to AA` = "YRI to AD",
            `AA to YRI` = "AD to YRI",
            `CEU to AA` = "CEU to AD",
            `AA to CEU` = "AD to CEU",
            `YRI to YRI` = "YRI to YRI",
            `CEU to CEU` = "CEU to CEU"
        )
    ) %>%
    group_by(Train_Test, Seed, Original_Model, YRI_proportion) %>%
    summarize(
        Post_Hoc_Power = sum(P_value < 0.05/ngenes),
        h2 = mean(Heritability)
    )
x.power2 = x.power2 %>%
    group_by(Original_Model) %>%
    summarize(Heritability = paste0("Heritability = ", round(mean(h2), 2))) %>%
    merge(., x.power2, by = "Original_Model")

# this applies a light grey color to pop-into-itself scenarios
my.colors.withgrey.admix = c(
    "CEU to AD"  = "blue",
    "YRI to AD"  = "red",
    "AD to CEU"  = "black",
    "AD to YRI"  = "black",
    "YRI to CEU" = "lightgrey",
    "CEU to YRI" = "lightgrey"
)
my.linetypes.admix = c(
    "CEU to AD"  = "solid",
    "YRI to AD"  = "solid",
    "AD to CEU"  = "dashed",
    "AD to YRI"  = "dotted",
    "YRI to CEU" = "dashed",
    "CEU to YRI" = "dotted"
)
my.colors.admix = c(
    "CEU to AD"  = "blue",
    "YRI to AD"  = "red",
    "AD to CEU"  = "black",
    "AD to YRI"  = "black",
    "YRI to CEU" = "red",
    "CEU to YRI" = "blue"
)
my.shapes.admix = c(
    "CEU to AD"  = "circle",
    "YRI to AD"  = "circle",
    "AD to CEU"  = "square",
    "AD to YRI"  = "triangle",
    "YRI to CEU" = "square",
    "CEU to YRI" = "triangle"
)
my.breaks.diffpop.admix = c(
    "YRI to AD",
    "AD to YRI",
    "AD to CEU",
    "CEU to AD",
    "YRI to CEU",
    "CEU to YRI"
)
x.power.admix.summary = x.power2  %>%
    group_by(Train_Test, Original_Model, YRI_proportion) %>%
    summarize(
        Power = sum(Post_Hoc_Power) / nseeds,
        Power_SE = sd(Post_Hoc_Power) / sqrt(n()),
        Power_CI_lo = Power - 1.96*Power_SE,
        Power_CI_hi = Power + 1.96*Power_SE,
        h2_Mean = mean(h2),
        h2_SE = sd(h2) / sqrt(n())
    ) %>%
    as.data.table

cat("plotting g9\n\n")
g9 = ggplot(x.power2, aes(x = YRI_proportion, y = Post_Hoc_Power, color = Train_Test)) +
        geom_smooth(aes(linetype = Train_Test), method = "glm", method.args = list(family = "binomial"), se = TRUE, alpha = 0.1, span = 0.5) +
        scale_shape_manual(
            name   = "Train to Test",
            values = my.shapes.admix,
            breaks = my.breaks.diffpop.admix
        ) +
        scale_linetype_manual(
            name   = "Train to Test",
            values = my.linetypes.admix,
            breaks = my.breaks.diffpop.admix
        )  +
        scale_color_manual(
            name   = "Train to Test",
            values = my.colors.withgrey.admix,
            breaks = my.breaks.diffpop.admix
        ) +
        facet_grid(~ Heritability) +
        theme_klk() +
        theme(legend.key.width = unit(1, "in")) + 
        guides(
            fill     = guide_legend(keywidth = 1, keyheight = 1),
            linetype = guide_legend(keywidth = 6, keyheight = 1),
            color    = guide_legend(keywidth = 6, keyheight = 1)
        ) +
        ylab("Power") +
        xlab("YRI proportion") +
        ggtitle(
            bquote("Power of TWAS association tests for varying proportions of YRI admixture"),
            subtitle = bquote(k~" = 10 eQTLs per gene, 50% shared eQTLs between populations") 
        )
g9.path = file.path(output.dir, paste("twas.sim.results.power.curves.admixvary", plot.type, sep = "."))
ggsave(g9, file = g9.path, units = "in", width = 18, height = 6)

# make supplementary tables with raw power estimates for varying admixture
admix.table.path = file.path(output.dir, "admixvary.power.supptables.txt") 
x.power.admix.summary %>%
    dplyr::filter(Original_Model %in% admix.effect.sizes) %>%
    mutate(Power_CI = paste(round(Power_CI_lo, 3), round(Power_CI_hi, 3), sep = " - ")) %>%
    dplyr::select(Train_Test, Original_Model, YRI_proportion, Power, Power_CI) %>%
    as.data.table %>%
    arrange(., Original_Model, Train_Test, YRI_proportion) %>%
    fwrite(., file = admix.table.path, sep = "\t", quote = FALSE)
