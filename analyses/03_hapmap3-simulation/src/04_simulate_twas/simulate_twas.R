suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(assertthat))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))

# parse command line variables
option_list = list(
    make_option(
        c("-a", "--results-directory"),
        type    = "character",
        default = NULL, 
        help    = "The directory path to the results stored in *Rdata files.", 
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
        c("-c", "--output-file-prefix"),
        type    = "character",
        default = NULL, 
        help    = "File name for results file, saved as a data.table.",
        metavar = "character"
    ),
    make_option(
        c("-d", "--num-eQTLs"),
        type    = "integer",
        default = NULL, 
        help    = "File name for results file, saved as a data.table.",
        metavar = "character"
    ),
    make_option(
        c("-e", "--proportion-shared-eQTL"),
        type    = "numeric",
        default = NULL, 
        help    = "The proportion of eQTLs shared between populations.",
        metavar = "numeric"
    ),
    make_option(
        c("-f", "--same-eQTLs"),
        type    = "logical",
        default = NULL, 
        help    = "Do populations share all eQTLs?", 
        metavar = "numeric"
    ),
    make_option(
        c("-g", "--same-eQTL-effects"),
        type    = "logical",
        default = NULL, 
        help    = "Do populations share eQTL effect sizes?",
        metavar = "logical"
    ),
    make_option(
        c("-i", "--random-seed"),
        type    = "integer",
        default = "2018", 
        help    = "Random seed for reproducible simulatino results. [default = %default]", 
        metavar = "integer"
    ),
#    make_option(
#        c("-j", "--h2-TWAS"),
#        type    = "numeric",
#        default = 0.3, 
#        help    = "The heritability of the simulated TWAS phenotype. [default = %default]", 
#        metavar = "numeric"
#    ),
    make_option(
        c("-k", "--same-TWAS-genes"),
        type    = "logical",
        default = TRUE, 
        help    = "Do populations share the same causal TWAS genes? [default = %default]", 
        metavar = "logical"
    ),
    make_option(
        c("-l", "--same-TWAS-effects"),
        type    = "logical",
        default = TRUE, 
        help    = "Do populations share the same effect sizes at causal TWAS genes? [default = %default]", 
        metavar = "logical"
    ),
    make_option(
        c("-m", "--TWAS-noise-mean"),
        type    = "numeric",
        default = 0.0, 
        help    = "The mean of environmental noise of the TWAS phenotype. [default = %default]", 
        metavar = "numeric"
    ),
    make_option(
        c("-n", "--TWAS-noise-standard-deviation"),
        type    = "numeric",
        default = 0.1, 
        help    = "The standard deviation of the environmental variance of the TWAS phenotype. [default = %default]", 
        metavar = "logical"
    ),
    make_option(
        c("-o", "--num-genes"),
        type    = "integer",
        default = 98, 
        help    = "The total number of genes available in the TWAS. [default = %default]", 
        metavar = "integer"
    ),
    make_option(
        c("-p", "--num-samples"),
        type    = "integer",
        default = 1000, 
        help    = "The number of samples in each population (should be the same for all pops).[default = %default]", 
        metavar = "integer"
    ),
    make_option(
        c("-q", "--num-causal-genes"),
        type    = "integer",
        default = 1, 
        help    = "The number of causal genes in the TWAS.[default = %default]", 
        metavar = "integer"
    ),
    make_option(
        c("-r", "--plot-type"),
        type    = "character",
        default = "png", 
        help    = "The plot extension (e.g. 'png', 'pdf', etc.) for output plots .[default = %default]", 
        metavar = "character"
    ),
    make_option(
        c("-s", "--effect-size"),
        type    = "numeric",
        default = "0.0", 
        help    = "The effect size for the causal genes .[default = %default]", 
        metavar = "numeric"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

rdata.file.dir    = opt$results_directory 
output.dir        = opt$output_directory
output.filepfx    = opt$output_file_prefix

num.eqtls         = as.integer(opt$num_eQTLs)
prop.shared.eqtl  = as.numeric(opt$proportion_shared_eQTL)
same.eqtls        = as.logical(opt$same_eQTLs) 
same.effects      = as.logical(opt$same_eQTL_effects) 

seed              = as.numeric(opt$random_seed) 
#h2.twas           = as.numeric(opt$h2_TWAS) 
same.twas.genes   = as.logical(opt$same_TWAS_genes)
same.twas.effects = as.logical(opt$same_TWAS_effects)
ngenes            = as.integer(opt$num_genes)
nsamples          = as.integer(opt$num_samples)
ncausal.genes     = as.integer(opt$num_causal_genes)
plot.type         = opt$plot_type 
effect.size       = as.numeric(opt$effect_size)
twas.noise.mean   = as.numeric(opt$TWAS_noise_mean)
twas.noise.sd     = as.numeric(opt$TWAS_noise_standard_deviation)



# ==========================================================================================
# script variables
# ==========================================================================================


twas.output.path      = file.path(output.dir, paste(output.filepfx, "txt", sep = ".")) 
twas.plot.output.path = file.path(output.dir, paste(output.filepfx, plot.type, sep = ".")) 


# load current simulated gene expression dataset
#load("./ELFN2_simulation_prediction_admixedpop_sameeQTLsFALSE_sameeffectsTRUE_k10_seed2018_propsharedeQTLs0.1.Rdata")
my.pattern = paste0("*simulation_prediction_admixedpop_sameeQTLs", same.eqtls, "_sameeffects", same.effects, "_k", num.eqtls, "_seed", seed, "_propsharedeQTLs", prop.shared.eqtl, ".Rdata")


rdata.files = list.files(rdata.file.dir, pattern = my.pattern, full.names = TRUE)
ndatafiles = length(rdata.files)
assert_that(ndatafiles == ngenes, msg = paste0("Number of genes ", ngenes, " does not match number of Rdata files ", ndatafiles, "\n")) 

twas.simulation = vector(mode = "list", length = ndatafiles)

for (my.file in rdata.files) {

    # load results for this configuration
    cat("Reading file ", my.file, "\n")
    load(my.file)

    # add original simulated expression values 
    twas.simulation$original.expression.1[[gene]]      = results.this.k$phenotypes$pop1 
    twas.simulation$original.expression.2[[gene]]      = results.this.k$phenotypes$pop2 
    twas.simulation$original.expression.admix[[gene]]  = results.this.k$phenotypes$admixpop

    # add predicted simulated expression values
    # must construct these from genotypes + predicted eQTL models
    twas.simulation$predicted.expression.1.to.1[[gene]]     = as.matrix(genos.1 %*% beta.pred.1)
    twas.simulation$predicted.expression.1.to.2[[gene]]     = as.matrix(genos.2 %*% beta.pred.1)
    twas.simulation$predicted.expression.1.to.admix[[gene]] = as.matrix(genos.admix %*% beta.pred.1)

    twas.simulation$predicted.expression.2.to.1[[gene]]     = as.matrix(genos.1 %*% beta.pred.2)
    twas.simulation$predicted.expression.2.to.2[[gene]]     = as.matrix(genos.2 %*% beta.pred.2)
    twas.simulation$predicted.expression.2.to.admix[[gene]] = as.matrix(genos.admix %*% beta.pred.2)

    twas.simulation$predicted.expression.admix.to.1[[gene]]     = as.matrix(genos.1 %*% beta.pred.admix)
    twas.simulation$predicted.expression.admix.to.2[[gene]]     = as.matrix(genos.2 %*% beta.pred.admix)
    twas.simulation$predicted.expression.admix.to.admix[[gene]] = as.matrix(genos.admix %*% beta.pred.admix)

}

# compile lists of expression vectors into matrices
original.expression.1      = as.data.frame(twas.simulation$original.expression.1)
original.expression.2      = as.data.frame(twas.simulation$original.expression.2)
original.expression.admix  = as.data.frame(twas.simulation$original.expression.admix)

predicted.expression.1.to.1     = as.data.frame(twas.simulation$predicted.expression.1.to.1)
predicted.expression.1.to.2     = as.data.frame(twas.simulation$predicted.expression.1.to.2)
predicted.expression.1.to.admix = as.data.frame(twas.simulation$predicted.expression.1.to.admix)

predicted.expression.2.to.1     = as.data.frame(twas.simulation$predicted.expression.2.to.1)
predicted.expression.2.to.2     = as.data.frame(twas.simulation$predicted.expression.2.to.2)
predicted.expression.2.to.admix = as.data.frame(twas.simulation$predicted.expression.2.to.admix)

predicted.expression.admix.to.1 = as.data.frame(twas.simulation$predicted.expression.admix.to.1)
predicted.expression.admix.to.2 = as.data.frame(twas.simulation$predicted.expression.admix.to.2)
predicted.expression.admix.to.admix = as.data.frame(twas.simulation$predicted.expression.admix.to.admix)

# set up simulated TWAS models 
twas.model.1   = sample.int(ngenes, size = ncausal.genes, replace = FALSE)
beta.twas.1    = matrix(0, ngenes, 1)
twas.effects.1 = effect.size 
beta.twas.1[twas.model.1] = twas.effects.1

# specify other two pops based on same.twas.genes and same.twas.effects
beta.twas.2     = matrix(0, ngenes, 1)
beta.twas.admix = matrix(0, ngenes, 1)

## should pop2 and admixed pop use same eQTL effect sizes?
## if not, then simulate new ones
#if ( same.twas.effects ) {
#    twas.effects.2     = twas.effects.1
#    twas.effects.admix = twas.effects.1
#} else { 
    twas.effects.2     = effect.size 
    twas.effects.admix = effect.size 
#}

# populations should share causal genes
twas.model.2     = twas.model.1
twas.model.admix = twas.model.1

# specify models for pop 2 and admixed pop
beta.twas.2[twas.model.2]         = twas.effects.2
beta.twas.admix[twas.model.admix] = twas.effects.admix

# specify the phenotype/environmental noise
# this is parametrized by the desired heritability h2 and # of eQTLs k
# only works if 0 <= h2 <= 1 !!! 
twas.pheno.1     = as.matrix(original.expression.1) %*% beta.twas.1
twas.pheno.2     = as.matrix(original.expression.2) %*% beta.twas.2
twas.pheno.admix = as.matrix(original.expression.admix) %*% beta.twas.admix

#twas.pheno.sd.1     = sqrt( var( twas.pheno.1 ) * (1 - h2.twas) / h2.twas )
#twas.pheno.sd.2     = sqrt( var( twas.pheno.2 ) * (1 - h2.twas) / h2.twas )
#twas.pheno.sd.admix = sqrt( var( twas.pheno.admix ) * (1 - h2.twas) / h2.twas )

twas.pheno.var.1     = var(twas.pheno.1)
twas.pheno.var.2     = var(twas.pheno.2)
twas.pheno.var.admix = var(twas.pheno.admix)
pheno.noise          = matrix(rnorm(nsamples, twas.noise.mean, twas.noise.sd))

## using original simulated expression values, simulate a phenotype for each pop
#twas.y.1     = twas.pheno.1 + as.matrix(rnorm(nsamples, 0, twas.pheno.sd.1))
#twas.y.2     = twas.pheno.2 + as.matrix(rnorm(nsamples, 0, twas.pheno.sd.2))
#twas.y.admix = twas.pheno.admix + as.matrix(rnorm(nsamples, 0, twas.pheno.sd.admix))

twas.y.1     = twas.pheno.1 + pheno.noise 
twas.y.2     = twas.pheno.2 + pheno.noise 
twas.y.admix = twas.pheno.admix + pheno.noise 

twas.y.sd.1     = sd(twas.y.1)
twas.y.sd.2     = sd(twas.y.2)
twas.y.sd.admix = sd(twas.y.admix)

# subroutine to perform phenotype-expression regressions in one go 
perform.twas.regressions = function(phenotype, expression, from.pop, to.pop, original.model, stddev, h2.g, h2.e) {

    ### two subroutines for use with apply
    # this one performs the marginal regressions in one go
    perform.regression = function(z) {
        my.lm         = lm(phenotype ~ z)
        my.lm.summary = summary(my.lm)
        p.value       = my.lm.summary$coefficients[2,4]
        twas.effect   = my.lm.summary$coefficients[2,1]
        t.value       = my.lm.summary$coefficients[2,3]
        stderr        = my.lm.summary$coefficients[2,2]
        return(list("p.value" = p.value, "twas.effect" = twas.effect, "t.value" = t.value, "stderr" = stderr))
    }

    # this one computes the power for each effect size in one go
    compute.power = function(z, s) {
        n = dim(expression)[1]
        p = dim(expression)[2]
        r = power.t.test(power = NULL, alternative="two.sided", delta = z, n = n, sd = s, type="one.sample", sig.level=0.05/p)$power
        return(r)
    }

    # perform regressions 
    my.regressions = apply(expression, 2, perform.regression)

    # compute power for each gene (both causal and noncausal) 
    my.power = apply(original.model, 1, function(z) compute.power(z, stddev))
    
    # save gene names to link to regression output
    my.genes = colnames(expression)

    # reformat/rename the output 
    my.output = my.regressions %>%
        map(as.data.frame) %>%
        bind_rows() %>%
        mutate(
            "genes"               = my.genes,
            "from.pop"            = from.pop,
            "to.pop"              = to.pop,
            "seed"                = seed,
            "ncausal.genes"       = ncausal.genes,
            "same.causal.genes"   = same.twas.genes,
            "same.causal.effects" = same.twas.effects,
            "prop.shared.eqtl"    = prop.shared.eqtl,
            "num.causal.eqtl"     = k,
            "same.causal.eqtl"    = same.eqtls,
            "same.eqtl.effects"   = same.effects,
            "original.model"      = original.model,
            "power.to.detect"     = my.power,
            "pheno.var"           = stddev^2,
            "h2.genetic"          = h2.g,
            "h2.environmental"    = h2.e,
            "h2.overall"          = h2.g / (h2.g + h2.e)
        ) %>%
        select(genes, twas.effect, stderr, t.value, p.value, from.pop, to.pop, seed, ncausal.genes, same.causal.genes,
            same.causal.effects, prop.shared.eqtl, num.causal.eqtl, same.causal.eqtl, same.eqtl.effects, original.model,
            power.to.detect, pheno.var, h2.genetic, h2.environmental, h2.overall
        ) %>%
        rename("Gene_Name" = genes, "TWAS_Effect" = twas.effect, "StdErr" = stderr, "T_value" = t.value, "P_value" = p.value,
            "Train_Pop" = from.pop, "Test_Pop" = to.pop, "Seed" = seed, "Num_Causal_Genes" = ncausal.genes,
            "Same_Causal_Genes" = same.causal.genes, "Same_Causal_Effects" = same.causal.effects,
            "Prop_Shared_eQTL" = prop.shared.eqtl, "Num_Causal_eQTL" = num.causal.eqtl, "Same_Causal_eQTL" = same.causal.eqtl,
            "Same_eQTL_Effects" = same.eqtl.effects, "Original_Model" = original.model, "Power" = power.to.detect,
            "Phenotype_Variance" = pheno.var, "Genetic_Variance" = h2.genetic, "Environmental_Variance" = h2.environmental,
            "Heritability" = h2.overall
        )
    return(data.table(my.output))
} 

# calculate a variance once to pass to regression functions
twas.noise.var = twas.noise.sd^2

# simulate TWAS regressions!
twas.1.to.1     = perform.twas.regressions(twas.y.1, predicted.expression.1.to.1, "CEU", "CEU", beta.twas.1, twas.y.sd.1, twas.pheno.var.1, twas.noise.var) 
twas.2.to.1     = perform.twas.regressions(twas.y.1, predicted.expression.2.to.1, "YRI", "CEU", beta.twas.1, twas.y.sd.1, twas.pheno.var.1, twas.noise.var) 
twas.admix.to.1 = perform.twas.regressions(twas.y.1, predicted.expression.admix.to.1, "AA", "CEU", beta.twas.1, twas.y.sd.1, twas.pheno.var.1, twas.noise.var) 

twas.1.to.2     = perform.twas.regressions(twas.y.2, predicted.expression.1.to.2, "CEU", "YRI", beta.twas.2, twas.y.sd.2, twas.pheno.var.2, twas.noise.var) 
twas.2.to.2     = perform.twas.regressions(twas.y.2, predicted.expression.2.to.2, "YRI", "YRI", beta.twas.2, twas.y.sd.2, twas.pheno.var.2, twas.noise.var) 
twas.admix.to.2 = perform.twas.regressions(twas.y.2, predicted.expression.admix.to.2, "AA", "YRI", beta.twas.2, twas.y.sd.2, twas.pheno.var.2, twas.noise.var) 

twas.1.to.admix     = perform.twas.regressions(twas.y.admix, predicted.expression.1.to.admix, "CEU", "AA", beta.twas.admix, twas.y.sd.admix, twas.pheno.var.admix, twas.noise.var) 
twas.2.to.admix     = perform.twas.regressions(twas.y.admix, predicted.expression.2.to.admix, "YRI", "AA", beta.twas.admix, twas.y.sd.admix, twas.pheno.var.admix, twas.noise.var) 
twas.admix.to.admix = perform.twas.regressions(twas.y.admix, predicted.expression.admix.to.admix, "AA", "AA", beta.twas.admix, twas.y.sd.admix, twas.pheno.var.admix, twas.noise.var) 

twas.results = rbindlist(list(twas.1.to.1, twas.2.to.1, twas.admix.to.1, twas.1.to.2, twas.2.to.2, twas.admix.to.2,
    twas.1.to.admix, twas.2.to.admix, twas.admix.to.admix))

# we only care about causal effects
# purge noncausal genes and write result to file
twas.results = twas.results %>% dplyr::filter(Original_Model != 0)
fwrite(x = twas.results, file = twas.output.path, sep = "\t", quote = FALSE, na = "NA")

# save output
save.image(file.path(output.dir, paste(output.filepfx, "Rdata", sep = ".")))
