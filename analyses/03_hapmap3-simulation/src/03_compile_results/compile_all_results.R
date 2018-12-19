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
        c("-r", "--results-directory"),
        type    = "character",
        default = NULL, 
        help    = "The directory path to the results stored in *Rdata files.", 
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
        default = "1kg.all.results.txt", 
        help    = "File name for results file, saved as a data.table. [default = %default]", 
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
rdata.file.dir  = opt$results_directory 
output.dir      = opt$output_directory
output.filename = opt$output_filename

output.txt.path = file.path(output.dir, output.filename) 


# ==========================================================================================
# script code 
# ==========================================================================================

# load all results files
rdata.files = list.files(rdata.file.dir, pattern = "*.Rdata", full.names = TRUE)

# how many Rdata files will we parse?
cat("Parsing ", length(rdata.files), " results files...\n")

# preallocate list for storing data frames
all.results.list = vector(mode = "list", length = length(rdata.files)) 

# keep counter for indexing
i = 0
for (my.file in rdata.files) {

    # increment counter
    i = i + 1

    # load results for this configuration
    cat("Reading file ", my.file, "\n")
    load(my.file)

    # make a molten data frame from R2 results=
    # will repeat for correlations and merge data tables
    my.r2.results   = melt(results.this.k$r2, id.vars = "test_pop", variable.name = "train_pop", value.name = "R2")
    my.corr.results = melt(results.this.k$corr, id.vars = "test_pop", variable.name = "train_pop", value.name = "Correlation") 

#    # compile information theoretic results
#    # these require some care; will only stick them to pop --> pop prediction rows, since they don't make sense otherwise
#    my.precision.results = melt(results.this.k$precision, value.name = "Precision") 
#    my.sensitivity.results = melt(results.this.k$sensitivity, value.name = "Sensitivity") 
#    my.specificity.results = melt(results.this.k$specificity, value.name = "Specificity") 
#
#    colnames(my.precision.results)[2] = "train_pop"
#    colnames(my.sensitivity.results)[2] = "train_pop"
#    colnames(my.specificity.results)[2] = "train_pop"
#
#    my.precision.results$test_pop = my.precision.results$train_pop
#    my.sensitivity.results$test_pop = my.precision.results$train_pop
#    my.specificity.results$test_pop = my.precision.results$train_pop

    my.results = data.table(merge(my.r2.results, my.corr.results, by = c("test_pop", "train_pop"), all = TRUE))
    #my.results = list(my.r2.results, my.corr.results, my.precision.results, my.sensitivity.results, my.specificity.results) %>%
#        reduce(full_join, by = c("train_pop", "test_pop"))

    # add current scenario information. there are four scenario variables:
    # same.eQTLs = T or F
    # same.eQTL.effects = T or F
    # k = 1, 5, 10, or 20
    # frac.same.eQTLs = 0.0, 0.1, 0.2, ... 
    # seed is the random seed
    my.results$k = k
    my.results$same_eqtls = same.eQTLs
    my.results$same_effects = same.eQTL.betas
    my.results$gene = gene
    my.results$prop_shared_eqtl = frac.same.eQTLs
    my.results$seed = seed

    # append results to list
    all.results.list[[i]] = my.results
}

# combine all data.tables in one go
all.results = rbindlist(all.results.list, use.names = TRUE)

# want to rename populations:
# pop1 --> CEU
# pop2 --> YRI
# pop3 --> AA
# can tinker with levels of factor variables to get order correct
all.results$Train_Pop = as.character(all.results$train_pop)
all.results$Train_Pop[all.results$train_pop == "pop1"] = "CEU"
all.results$Train_Pop[all.results$train_pop == "pop2"] = "YRI"
all.results$Train_Pop[all.results$train_pop == "admixpop"] = "AA"
all.results$Train_Pop = as.factor(all.results$Train_Pop)

all.results$Test_Pop = as.character(all.results$test_pop)
all.results$Test_Pop[all.results$test_pop == "pop1"] = "CEU"
all.results$Test_Pop[all.results$test_pop == "pop2"] = "YRI"
all.results$Test_Pop[all.results$test_pop == "admixpop"] = "AA"
all.results$Test_Pop = as.factor(all.results$Test_Pop)

# save compiled results   
fwrite(all.results, file = output.txt.path, quote = FALSE, sep = "\t") 
