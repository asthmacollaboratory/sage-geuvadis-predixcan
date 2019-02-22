#!/usr/bin/env Rscript --vanilla
# want to see how genetic distance varies with prediction accuracy in GEUVADIS
# simple plot of F_ST versus R2

# load libraries
library(data.table)
library(ggplot2)

# enable plotting for PNG
options(bitmapType='cairo')

# load data
x = fread("../../datafiles/geuvadis.fst.r2.txt")

# shouldn't need this unless we want to jitter the points
# no harm in setting it anyways
set.seed(2019)

# construct simple point plot
g = ggplot(x, aes(x = R2, y = Fst, color = Imputation_Direction, shape = Imputation_Direction)) +
    geom_point(size = 3) + ggtitle(bquote("Fixation index" ~ F[st] ~ "versus prediction accuracy" ~ R^{2})) +
    xlab(bquote(R^{2})) +
    ylab(bquote(F[st])) +
    theme(legend.title = element_blank())

# save to analysis folder
ggsave(filename = "../../analysis/geuvadis.fst.vs.r2.png", plot = g, dpi = 300)
