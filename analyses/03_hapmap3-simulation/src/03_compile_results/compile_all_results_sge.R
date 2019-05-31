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
suppressMessages(library(purrr))
suppressMessages(library(methods))
suppressMessages(library(optparse))


# parse command line variables
option_list = list(
    make_option(
        c("-r", "--results-file"),
        type    = "character",
        default = NULL,
        help    = "The file path to one *Rdata file.",
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
        c("-f", "--output-filename"),
        type    = "character",
        default = NULL,
        help    = "File name for results file, saved as a data.table.",
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
my.file         = opt$results_file
output.dir      = opt$output_directory
output.filename = opt$output_filename

output.txt.path = file.path(output.dir, output.filename)


# ==========================================================================================
# script code
# ==========================================================================================

# load results for this configuration
cat("Reading file ", my.file, "\n")
load(my.file)

# make a molten data frame from R2 results
# will repeat for correlations and merge data tables
my.r2.results   = melt(results.this.k$r2, id.vars = "test_pop", variable.name = "train_pop", value.name = "R2")
my.corr.results = melt(results.this.k$corr, id.vars = "test_pop", variable.name = "train_pop", value.name = "Correlation")

my.results = data.table(merge(my.r2.results, my.corr.results, by = c("test_pop", "train_pop"), all = TRUE))

# add current scenario information. there are multiple scenario variables:
# same.eQTLs = T or F
# same.eQTL.effects = T or F
# k = 1, 5, 10, or 20
# frac.same.eQTLs = 0.0, 0.1, 0.2, ...
# seed is the random seed
# CEU_prop and YRI_prop are respective admixture proportions from CEU and YRI
my.results$k                = results.this.k$k
my.results$same_eqtls       = results.this.k$same.eqtls
my.results$same_effects     = results.this.k$same.effects
my.results$gene             = results.this.k$gene
my.results$prop_shared_eqtl = results.this.k$prop.shared.eqtl
my.results$seed             = results.this.k$seed
my.results$CEU_prop         = results.this.k$admix.pop1
my.results$YRI_prop         = results.this.k$admix.pop2

# want to rename populations:
# pop1 --> CEU
# pop2 --> YRI
# pop3 --> AA
# can tinker with levels of factor variables to get order correct
my.results$Train_Pop = as.character(my.results$train_pop)
my.results$Train_Pop[my.results$train_pop == "pop1"] = "CEU"
my.results$Train_Pop[my.results$train_pop == "pop2"] = "YRI"
my.results$Train_Pop[my.results$train_pop == "admixpop"] = "AA"
my.results$Train_Pop = as.factor(my.results$Train_Pop)

my.results$Test_Pop = as.character(my.results$test_pop)
my.results$Test_Pop[my.results$test_pop == "pop1"] = "CEU"
my.results$Test_Pop[my.results$test_pop == "pop2"] = "YRI"
my.results$Test_Pop[my.results$test_pop == "admixpop"] = "AA"
my.results$Test_Pop = as.factor(my.results$Test_Pop)

# save compiled results
fwrite(my.results, file = output.txt.path, quote = FALSE, sep = "\t")
