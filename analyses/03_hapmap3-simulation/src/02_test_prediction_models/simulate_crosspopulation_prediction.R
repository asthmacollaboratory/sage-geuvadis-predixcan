###!/usr/bin/env Rscript --vanilla
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
# This script simulates phenotypes and computes predictive models in three populations.
# The populations are specified by genotype files from which phenotypes are simulated.
# Each phenotype contains a prespecified number of true eQTLs; other genotypes have 0 effect.
# Each phenotype/genotype set is then used to train a predictive model.
# The training scheme uses elastic net regression with a nested crossvalidation scheme.
# ==========================================================================================

# ==========================================================================================
# script options
# ==========================================================================================

# don't emit warnings, only errors
options(warn = -1)

# ==========================================================================================
# libraries
# ==========================================================================================
suppressMessages(library(glmnet))
suppressMessages(library(methods))
suppressMessages(library(assertthat))
suppressMessages(library(data.table))
suppressMessages(library(doParallel))
suppressMessages(library(optparse))


# parse command line variables
option_list = list(
    make_option(
        c("-g", "--gene-name"),
        type    = "character",
        default = NULL,
        help    = "The name of the gene being analyzed, used for prefixing output",
        metavar = "character"
    ),
    make_option(
        c("-g1", "--genotypes-pop1"),
        type    = "character",
        default = NULL,
        help    = "PLINK RAW file for population 1",
        metavar = "character"
    ),
    make_option(
        c("-g2", "--genotypes-pop2"),
        type    = "character",
        default = NULL,
        help    = "PLINK RAW file for population 2",
        metavar = "character"
    ),
    make_option(
        c("-ga", "--genotypes-admix"),
        type    = "character",
        default = NULL,
        help    = "PLINK RAW file for admixed population",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-directory"),
        type    = "character",
        default = NULL,
        help    = "Directory where prediction results will be stored in Rdata format",
        metavar = "character"
    ),
    make_option(
        c("-seq", "--same-eqtls"),
        type    = "logical",
        default = TRUE,
        help    = "should predictive models use exact same eQTL positions? [default = %default]",
        metavar = "logical"
    ),
    make_option(
        c("-f", "--fraction-overlapping-eqtls"),
        type    = "double",
        default = 0.0,
        help    = "what fraction of eQTL positions should pop1 and pop2 share? Only applies when --same-eqtls=FALSE [default = %default]",
        metavar = "double"
    ),
    make_option(
        c("-sef", "--same-eqtl-effects"),
        type    = "logical",
        default = TRUE,
        help    = "should predictive models use same eQTL effects? [default = %default]",
        metavar = "logical"
    ),
    make_option(
        c("-k", "--num-eqtls"),
        type    = "integer",
        default = "1",
        help    = "number of eQTLs to simulate [default= %default]",
        metavar = "integer"
    ),
    make_option(
        c("-ne", "--nfolds-external"),
        type    = "integer",
        default = "5",
        help    = "how many folds should be run for the external crossvalidation? [default= %default]",
        metavar = "integer"
    ),
    make_option(
        c("-ni", "--nfolds-internal"),
        type    = "integer",
        default = "10",
        help    = "how many folds should be run for the internal nested crossvalidation? [default= %default]",
        metavar = "integer"
    ),
    make_option(
        c("-np", "--nfolds-parallel"),
        type    = "integer",
        default = "4",
        help    = "how many folds should be run in parallel? [default= %default]",
        metavar = "integer"
    ),
    make_option(
        c("-rs", "--random-seed"),
        type    = "integer",
        default = "2018",
        help    = "random seed for reproducibility [default= %default]",
        metavar = "integer"
    ),
    make_option(
        c("-em", "--eqtl-mean"),
        type    = "double",
        default = "0.0",
        help    = "mean of eQTL effect sizes [default= %default]",
        metavar = "double"
    ),
    make_option(
        c("-es", "--eqtl-sd"),
        type    = "double",
        default = "1.0",
        help    = "standard deviation of eQTL effect sizes [default= %default]",
        metavar = "double"
    ),
    make_option(
        c("-h2", "--heritability"),
        type    = "double",
        default = "0.15",
        help    = "heritability of simulated phenotypes[default= %default]",
        metavar = "double"
    ),
    make_option(
        c("-a", "--alpha"),
        type    = "double",
        default = "0.5",
        help    = "alpha mixing parameter for elastic net [default= %default]",
        metavar = "double"
    ),
    make_option(
        c("-m", "--maf-min"),
        type    = "double",
        default = "0.05",
        help    = "minimum admissible MAF for genotypes used in model construction [default= %default]",
        metavar = "double"
    ),
    make_option(
        c("-M", "--maf-max"),
        type    = "double",
        default = "0.95",
        help    = "maximum admissible MAF for genotypes used in model construction [default= %default]",
        metavar = "double"
    ),
    make_option(
        c("-A1", "--admix-proportion-pop1"),
        type    = "double",
        default = "0.2",
        help    = "Proportion of haplotypes in admixed population drawn from population 1[default= %default]",
        metavar = "double"
    ),
    make_option(
        c("-A2", "--admix-proportion-pop2"),
        type    = "double",
        default = "0.8",
        help    = "Proportion of haplotypes in admixed population drawn from population 2[default= %default]",
        metavar = "double"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n\n")
print(opt)

cat("This script will look for libraries here:\n")
print(.libPaths())
cat("Missing functions (e.g. `get_init_params` in glmnet) suggest a failure to correctly load a library!\n\n")  ## <-- THIS ERROR IS EVIL


# ==========================================================================================
# script variables
# ==========================================================================================
gene             = opt$gene_name                   # name of the gene being analyzed
genos.1.path     = opt$genotypes_pop1              # filepath to PLINK RAW file for pop 1 (e.g. CEU)
genos.2.path     = opt$genotypes_pop2              # filepath to PLINK RAW file for pop 2 (e.g. YRI)
genos.admix.path = opt$genotypes_admix             # filepath to PLINK RAW file for pop 2 (e.g. AA)
output.dir       = opt$output_directory            # where will the prediction output be stored?
eqtl.mean        = opt$eqtl_mean                   # mean of eQTL effect sizes
eqtl.sd          = opt$eqtl_sd                     # sd of eQTL effect sizes
h2               = opt$heritability                # heritability of phenotype, which parametrizes phenotypic noise
alpha            = opt$alpha                       # elastic net mixing parameter
seed             = opt$random_seed                 # random seed for reproducibility
k                = as.numeric(opt$num_eqtls)       # number of true eQTLs to simulate
same.eQTLs       = opt$same_eqtls                  # should the predictive models use the exact same eQTLs?
frac.same.eQTLs  = opt$fraction_overlapping_eqtls  # what fraction of the eQTLs should the three populations share?
same.eQTL.betas  = opt$same_eqtl_effects           # should the predictive models use the same eQTL effect sizes?
nfolds.external  = opt$nfolds_external             # how many external crossvalidation folds should be run?
nfolds.internal  = opt$nfolds_internal             # how many internal (nested) CV folds should be run?
nfolds.parallel  = opt$nfolds_parallel             # how many folds should be run in parallel?
maf.min          = opt$maf_min                     # minimum admissible minor allele frequency from the genotype files
maf.max          = opt$maf_max                     # maximum admissible minor allele frequency from the genotype files
admix.prop.pop1  = as.numeric(opt$admix_proportion_pop1) # proportion of haplotypes from pop1 (e.g. CEU, 20%)
admix.prop.pop2  = as.numeric(opt$admix_proportion_pop2) # proportion of haplotypes from pop2 (e.g. YRI, 80%)

# don't forget to set the seed!
set.seed(seed)

# ==========================================================================================
# subroutines
# ==========================================================================================

# compute the minor allele frequency for one column of a genotype dosage matrix
# checks that all entries of the column are one of 0, 1, 2, or NA
# throws an error otherwise
maf.onesnp = function(x) {
    n   = length(x)
    n0  = sum(x == 0, na.rm = TRUE)
    n1  = sum(x == 1, na.rm = TRUE)
    n2  = sum(x == 2, na.rm = TRUE)
    nna = sum(is.na(x), na.rm = TRUE)
    assert_that( n0 + n1 + n2 + nna == n)
    return( sum(x, na.rm = TRUE) / (2*n - nna) )
}

# compute minor allele frequencies for all columns of a genotype dosage matrix
maf = function(x) {
    y = apply(x, 2, maf.onesnp)
    return(y)
}

# impute missing entries of a genotype matrix using precomputed minor allele frequencies
impute.with.maf = function(x, mafs) {
    z = sapply(1:dim(x)[2], function(i) {y = x[, ..i]; y[is.na(y)] = 2*mafs[i]; return(y)} )
    z = as.data.table(z)
    colnames(z) = colnames(x)
    return(z)
}

# subroutine to partition a range of numbers
# this yields indices to pull held-out-samples
chunk.samples = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

# train a predictive model for one population and one continuous phenotype
# use elastic net with selectable alpha parameter
# train with a nested crossvalidation scheme:
# -- outer nest is "nfolds.external"-CV
# -- inner nest is "nfolds.internal"-CV
train.predictive.model = function(x, y, nfolds = length(y) - 1, verbose = TRUE, alpha = 0.5, nlambda = 100, parallel = TRUE, lambda.min.ratio = 0.01){

	n   = length(y)
    n.x = dim(x)[1]
    p   = dim(x)[2]

    # error checking
	assert_that(n.x == n)    # ensure that x and y have the same number of rows
    assert_that(nfolds <= n) # ensure that number of folds does not exceed number of samples

    # preallocate matrices to store output
	best.betas = matrix(-Inf, p, n)
	y.pred     = matrix(-Inf, n, 1)

    # make crossvalidation sample
    folds = chunk.samples(1:n, nfolds.external)

	# perform nested crossvalidation scheme
    # outer loop is leave-one-out
    # inner loop uses nfolds crossvalidation
	for (held.out.sample in folds) {

		# subset training/testing data
		X = x[-held.out.sample,]
		Y = y[-held.out.sample,]

		# cross-validate!
		#if (verbose) cat(paste0("\tcrossvalidating sample ", held.out.sample, " at timestamp ", Sys.time(), "\n"))
		if (verbose) cat(paste0("\tcrossvalidating ", length(held.out.sample), " held-out-samples at timestamp ", Sys.time(), "\n"))

        my.glmnet = cv.glmnet(
            x = X,
            y = Y,
            nfolds   = nfolds.internal,
            family   = "gaussian",
            alpha    = alpha,
            keep     = TRUE,
            grouped  = FALSE,
            dfmax    = n - 1,
            pmax     = n - 1,
            lambda.min.ratio = lambda.min.ratio,
            nlambda  = nlambda,
            parallel = parallel
        )

		if (verbose) cat(paste0("\tdone at timestamp ", Sys.time(), "\n"))

		# parse results
        # find best lambda and corresponding model
		nrow.best    = which.min(my.glmnet$cvm) # position of best λ in cv.glmnet output
		my.best.beta = my.glmnet$glmnet.fit$beta[,nrow.best]
		best.betas[,held.out.sample] = my.best.beta

		# make prediction on held-out sample
		# can use s = lambda.min to use best predictive λ directly
		#y.heldout.expr.pred = predict(my.glmnet, newx = t(as.matrix(x[held.out.sample,])), s = "lambda.min")
		y.heldout.expr.pred = predict(my.glmnet, newx = as.matrix(x[held.out.sample,]), s = "lambda.min")

		y.pred[held.out.sample] = y.heldout.expr.pred
	}
	return(list("y" = y.pred, "betas" = best.betas))
}

load.geno.data = function(file.path){
    geno.data = fread(file.path, header = TRUE)[,-c(2:6)]  ## cols 2-5 not needed here
    colnames(geno.data)[1] = "SubjectID"
    setkey(geno.data, "SubjectID")
    setorderv(geno.data, "SubjectID", 1L)
    subject.ids = geno.data$SubjectID
    geno.data = geno.data[,-1]
    return(list("geno.data" = geno.data, "subject.ids" = subject.ids))
}

# ==========================================================================================
# load genotype data from file
# ==========================================================================================


cat(paste0("Start time: ", Sys.time(), "\n\n"))

cat("loading genotype data...\n")
pop.1       = load.geno.data(genos.1.path)
pop.2       = load.geno.data(genos.2.path)
pop.admix   = load.geno.data(genos.admix.path)
genos.1     = pop.1$geno.data
genos.2     = pop.2$geno.data
genos.admix = pop.admix$geno.data


# get sample sizes
# do this before culling genotypes by MAF
N.ancestral.1 = dim(genos.1)[1]
N.ancestral.2 = dim(genos.2)[1]
N.admix       = dim(genos.admix)[1]

# subsample the ancestral pops
# this power-matches their models
N.ancestral.min = min(N.ancestral.1, N.ancestral.2)
samples.1 = sample.int(N.ancestral.1, N.ancestral.min, replace = FALSE)
samples.2 = sample.int(N.ancestral.2, N.ancestral.min, replace = FALSE)

genos.1 = genos.1[samples.1, ]
genos.2 = genos.2[samples.2, ]


# get SNPs in common b/w all data
snps.1     = names(genos.1)
snps.2     = names(genos.2)
snps.admix = names(genos.admix)
commonsnps = intersect(snps.1, intersect(snps.2, snps.admix))

cat("Number of SNPs in common: ", length(commonsnps), "\n")

# subset genotype data to only SNPs present in all 3 pops
genos.1     = genos.1[, ..commonsnps]
genos.2     = genos.2[, ..commonsnps]
genos.admix = genos.admix[, ..commonsnps]


# compute MAFs of genotypes
# impute missing genotypes with 2*MAF

cat("Imputing missing dosages with 2*MAF...\n")
maf.1     = maf(genos.1)
maf.2     = maf(genos.2)
maf.admix = maf(genos.admix)

genos.1     = impute.with.maf(genos.1, maf.1)
genos.2     = impute.with.maf(genos.2, maf.2)
genos.admix = impute.with.maf(genos.admix, maf.admix)


# find SNPs within MAF limits
# then use them for subsetting
maf.pass.1     = maf.1 >= maf.min & maf.1 <= maf.max
maf.pass.2     = maf.2 >= maf.min & maf.2 <= maf.max
maf.pass.admix = maf.admix >= maf.min & maf.admix <= maf.max

maf.pass = maf.pass.1 & maf.pass.2 & maf.pass.admix

commonsnps.mafpass = commonsnps[maf.pass]

cat("Total number of common SNPs within MAF thresholds: ", length(commonsnps.mafpass), "\n")

genos.1     = as.matrix(genos.1[, ..commonsnps.mafpass])
genos.2     = as.matrix(genos.2[, ..commonsnps.mafpass])
genos.admix = as.matrix(genos.admix[, ..commonsnps.mafpass])


# get numbers of SNPs here. these numbers should be *exactly the same*
p.1     = dim(genos.1)[2]
p.2     = dim(genos.2)[2]
p.admix = dim(genos.admix)[2]
if ( !all(p.1 == p.2, p.2 == p.admix, p.admix == p.1) ) {
    stop("ERROR: Different SNPs between genotype files. Check SNPs in common between populations!\n\n")
}


# how many (common) SNPs are in this gene?
M = p.1
cat("M = ", M, "\n")
cat("p1 = ", p.1, "\n")
cat("p2 = ", p.2, "\n")
cat("p.admix = ", p.admix, "\n")

cat(N.ancestral.min, " samples used in each ancestral population.\n")
cat(M, " SNPs pass thresholds for analysis.\n")

# tidy the memory pool before simulations
gc()


# ==========================================================================================
# simulate phenotypes
# ==========================================================================================

## refix the random seed
#set.seed(seed)

# register a parallel backend
# this will enable glmnet to train models more quickly
registerDoParallel(nfolds.parallel)

cat("creating predictive models...\n")

cat(paste0("\nANALYSIS FOR NUMBER k = ", k, " eQTLs\n"))

snp.model.1         = sample.int(M, size = k, replace = FALSE)
beta.1              = matrix(0, M, 1)
eqtl.effects.1      = rnorm(k, eqtl.mean, eqtl.sd)
beta.1[snp.model.1] = eqtl.effects.1

# specify other two pops based on same.eQTLs and same.eQTL.betas
beta.2     = matrix(0, M, 1)
beta.admix = matrix(0, M, 1)

# should pop2 and admixed pop use same eQTL effect sizes?
# if not, then simulate new ones
if ( same.eQTL.betas ) {
    eqtl.effects.2     = eqtl.effects.1
    eqtl.effects.admix = eqtl.effects.1
} else {
    eqtl.effects.2     = rnorm(k, eqtl.mean, eqtl.sd)
    eqtl.effects.admix = rnorm(k, eqtl.mean, eqtl.sd)
}

# should pop2 and admixed pop use same eQTLs as pop 1?
# if not, then simulate new ones
if ( same.eQTLs ) {
    snp.model.2     = snp.model.1
    snp.model.admix = snp.model.1
} else {
    # slightly complicated sampling scheme here
    # simulate a fraction of eqtls in common b/w pop1, pop2
    # preserve those for pop2 and admix pop
    # for eqtls NOT in common, pop2 gets some random sample from elements of 1:M not in snp.model.1
    # this draws eqtls unique to pop2 from whatever eQTL positions are NOT in pop1
    # for admix pop, preserve eqtls in common b/w pop1 and pop2
    # then sample remaining eqtls from pop1, pop2 in ancestral fraction (20% of pop1, 80% of pop2)
    # note: the *ancestral* proportions of eqtls from each ancestral pop are fixed, but the % shared eqtls varies
    # if pop1, pop2 share no eqtls, then admix pop will still have eqtls in common with each ancestral pop
    eqtls.in.common = sample(snp.model.1, size = floor(frac.same.eQTLs * k), replace = FALSE)
    eqtls.different = sample(
        setdiff(1:M, snp.model.1),
        size    = ceiling((1 - frac.same.eQTLs) * k),
        replace = FALSE
    )
    snp.model.2     = c(eqtls.in.common, eqtls.different)
    snp.model.admix = c(eqtls.in.common,
        sample(snp.model.1, size = floor(admix.prop.pop1 * k), replace = FALSE),
        sample(snp.model.2, size = ceiling(admix.prop.pop2 * k), replace = FALSE)
    )
}

# specify models for pop 2 and admixed pop
beta.2[snp.model.2]         = eqtl.effects.2
beta.admix[snp.model.admix] = eqtl.effects.admix

# specify the phenotype/environmental noise
# this is parametrized by the desired heritability h2 and # of eQTLs k
# only works if 0 <= h2 <= 1 !!!
pheno.sd.1 = sqrt( var( genos.1 %*% beta.1 ) * (1 - h2) / h2 )
pheno.sd.2 = sqrt( var( genos.2 %*% beta.2 ) * (1 - h2) / h2 )
pheno.sd.admix = sqrt( var( genos.admix %*% beta.admix ) * (1 - h2) / h2 )

# using eQTL models, simulate a phenotype for each pop
y.1     = genos.1 %*% beta.1 + as.matrix(rnorm(N.ancestral.min, 0, pheno.sd.1))
y.2     = genos.2 %*% beta.2 + as.matrix(rnorm(N.ancestral.min, 0, pheno.sd.2))
y.admix = genos.admix %*% beta.admix + as.matrix(rnorm(N.admix, 0, pheno.sd.admix))

#print(y.1[1:5])
#print(y.2[1:5])
#print(y.admix[1:5])

# run elastic net on current configuration of y, beta
# run once for each population
cat("================= training model for pop 1 =================\n")
model.1 = train.predictive.model(genos.1, y.1, nfolds = nfolds.external, alpha = alpha, parallel = TRUE)

cat("================= training model for pop 2 =================\n")
model.2 = train.predictive.model(genos.2, y.2, nfolds = nfolds.external, alpha = alpha, parallel = TRUE)

cat("================= training model for admixed pop =================\n")
model.admix = train.predictive.model(genos.admix, y.admix, nfolds = nfolds.external, alpha = alpha, parallel = TRUE)

# now perform cross-population prediction
# consists of using prediction weights from one pop with genos from another pop
# average weights across validation sets, thereby collapsing into a single model from each pop
beta.pred.1     = apply(model.1[["betas"]], 1, mean, na.rm = TRUE)
beta.pred.2     = apply(model.2[["betas"]], 1, mean, na.rm = TRUE)
beta.pred.admix = apply(model.admix[["betas"]], 1, mean, na.rm = TRUE)

# predict with model from pop 1
y.2.from.1     = genos.2 %*% beta.pred.1
y.admix.from.1 = genos.admix %*% beta.pred.1

# compare quality of predictions
# use both R2 and correlations
r2.2.from.1       = summary(lm(y.2.from.1 ~ y.2))$r.squared
r2.admix.from.1   = summary(lm(y.admix.from.1 ~ y.admix))$r.squared
r2.1.from.1       = summary(lm(model.1[["y"]] ~ y.1))$r.squared

corr.2.from.1     = unname(cor.test(y.2.from.1, y.2, method = "spearman")$estimate)
corr.admix.from.1 = unname(cor.test(y.admix.from.1, y.admix, method = "spearman")$estimate)
corr.1.from.1     = unname(cor.test(model.1[["y"]], y.1, method = "spearman")$estimate)

# repeat for model from pop 2
y.1.from.2     = genos.1 %*% beta.pred.2
y.admix.from.2 = genos.admix %*% beta.pred.2

r2.1.from.2       = summary(lm(y.1.from.2 ~ y.1))$r.squared
r2.admix.from.2   = summary(lm(y.admix.from.2 ~ y.admix))$r.squared
r2.2.from.2       = summary(lm(model.2[["y"]] ~ y.2))$r.squared
corr.1.from.2     = unname(cor.test(y.1.from.2, y.1, method = "spearman")$estimate)
corr.admix.from.2 = unname(cor.test(y.admix.from.2, y.admix, method = "spearman")$estimate)
corr.2.from.2     = unname(cor.test(model.2[["y"]], y.2, method = "spearman")$estimate)

# repeat for model from admix pop
y.1.from.admix = genos.1 %*% beta.pred.admix
y.2.from.admix = genos.2 %*% beta.pred.admix

r2.1.from.admix       = summary(lm(y.1.from.admix ~ y.1))$r.squared
r2.2.from.admix       = summary(lm(y.2.from.admix ~ y.2))$r.squared
r2.admix.from.admix   = summary(lm(model.admix[["y"]] ~ y.admix))$r.squared
corr.1.from.admix     = unname(cor.test(y.1.from.admix, y.1, method = "spearman")$estimate)
corr.2.from.admix     = unname(cor.test(y.2.from.admix, y.2, method = "spearman")$estimate)
corr.admix.from.admix = unname(cor.test(model.admix[["y"]], y.admix, method = "spearman")$estimate)

# compile the results
r2.by.pop = data.table(
    "test_pop" = c("pop1", "pop2", "admixpop"),
    "pop1"     = c(r2.1.from.1, r2.2.from.1, r2.admix.from.1),
    "pop2"     = c(r2.1.from.2, r2.2.from.2, r2.admix.from.2),
    "admixpop" = c(r2.1.from.admix, r2.2.from.admix, r2.admix.from.admix)
)
corr.by.pop = data.table(
    "test_pop" = c("pop1", "pop2", "admixpop"),
    "pop1"     = c(corr.1.from.1, corr.2.from.1, corr.admix.from.1),
    "pop2"     = c(corr.1.from.2, corr.2.from.2, corr.admix.from.2),
    "admixpop" = c(corr.1.from.admix, corr.2.from.admix, corr.admix.from.admix)
)
predictive.models = list(
    "pop1"     = beta.pred.1,
    "pop2"     = beta.pred.2,
    "admixpop" = beta.pred.admix
)

# also compile the simulated models
original.models = list(
    "pop1"     = beta.1,
    "pop2"     = beta.2,
    "admixpop" = beta.admix
)
genotypes = list(
    "pop1"     = genos.1,
    "pop2"     = genos.2,
    "admixpop" = genos.admix
)
phenotypes = list(
    "pop1"     = y.1,
    "pop2"     = y.2,
    "admixpop" = y.admix
)

# compile information theoretic results
beta.1.nz = beta.1 != 0
beta.2.nz = beta.2 != 0
beta.admix.nz = beta.admix != 0

beta.1.zero = beta.1 == 0
beta.2.zero = beta.2 == 0
beta.admix.zero = beta.admix == 0

beta.pred.1.nz = beta.pred.1 != 0
beta.pred.2.nz = beta.pred.2 != 0
beta.pred.admix.nz = beta.pred.admix != 0

beta.pred.1.zero = beta.pred.1 == 0
beta.pred.2.zero = beta.pred.2 == 0
beta.pred.admix.zero = beta.pred.admix == 0

precision = list(
    "pop1"     = sum(beta.1.nz & beta.pred.1.nz) / sum(beta.pred.1.nz),
    "pop2"     = sum(beta.2.nz & beta.pred.2.nz) / sum(beta.pred.2.nz),
    "admixpop" = sum(beta.admix.nz & beta.pred.admix.nz) / sum(beta.pred.admix.nz)
)
sensitivity = list(
    "pop1"     = sum(beta.1.nz & beta.pred.1.nz) / k,
    "pop2"     = sum(beta.2.nz & beta.pred.2.nz) / k,
    "admixpop" = sum(beta.admix.nz & beta.pred.admix.nz) / k
)
specificity = list(
    "pop1"     = sum(beta.1.zero & beta.pred.1.zero) / sum(beta.1.zero),
    "pop2"     = sum(beta.2.zero & beta.pred.2.zero) / sum(beta.2.zero),
    "admixpop" = sum(beta.admix.zero & beta.pred.admix.zero) / sum(beta.admix.zero)
)

# compile list of all results and simulated models for current number of true eQTLs
results.this.k = list(
    "r2"   = r2.by.pop,
    "corr" = corr.by.pop,
    "k"    = k,
    "seed" = seed,
    "gene" = gene,
    "genotypes"    = genotypes,
    "phenotypes"   = phenotypes,
    "same.eqtls"   = same.eQTLs,
    "same.effects" = same.eQTL.betas,
    "precision"    = precision,
    "sensitivity"  = sensitivity,
    "specificity"  = specificity,
    "admix.pop1"   = admix.prop.pop1,
    "admix.pop2"   = admix.prop.pop2,
    "predictive.models" = predictive.models,
    "original.models"   = original.models,
    "prop.shared.eqtl"  = frac.same.eQTLs
)

cat("\nanalysis done.\n\n")


cat("saving results...\n")
datafile.path = file.path(output.dir, paste0(gene, "_simulation_prediction_admixedpop_sameeQTLs", same.eQTLs, "_sameeffects", same.eQTL.betas, "_k", k, "_propsharedeQTLs", frac.same.eQTLs, "_CEU", admix.prop.pop1, "_YRI", admix.prop.pop2, "_seed", seed, ".Rdata"))
save(results.this.k, file = datafile.path)

cat(paste0("End time: ", Sys.time(), "\n\n"))
