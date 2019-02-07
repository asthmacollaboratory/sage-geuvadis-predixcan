#!/usr/bin/env Rscript --vanilla
# ==========================================================================================
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys
#
# This script runs elastic net on the GEUVADIS RNA-Seq data in order to compute new PrediXcan weights.
#
# heavily modified from original PrediXcan script by Heather E. Wheeler (2015-02-02)
# https://github.com/hakyim/PrediXcan/blob/master/Paper-Scripts/Heather/DGN-calc-weights/01_imputedDGN-WB_CV_elasticNet.r
# in Github repo, see runscripts/run_01_imputedDGN-WB_CV_elasticNet_chr*sh and qsub.txt for tarbell job submission scripts
# ==========================================================================================

# ==========================================================================================
# load libraries
# ==========================================================================================

suppressMessages(library(methods))
suppressMessages(library(glmnet))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(caret))

# error tracking
options(show.error.locations = TRUE)

# main script function
# it's quite long and occupies most of the code in this file
# we will execute at end
compute.new.weights = function(){


    # parse command line variables
    option_list = list(
        make_option(
            c("-d", "--genotype-dosage-file"),
            type    = "character",
            default = NULL,
            help    = "File with genotype dosages for gene of interest.",
            metavar = "character"
        ),
        make_option(
            c("-e", "--expression-file"),
            type    = "character",
            default = NULL,
            help    = "File with gene expression measurements for gene of interest.",
            metavar = "character"
        ),
        make_option(
            c("-g", "--gene-name"),
            type    = "character",
            default = NULL,
            help    = "Name of gene being analyzed.",
            metavar = "character"
        ),
        make_option(
            c("-p", "--prediction-output"),
            type    = "character",
            default = NULL,
            help    = "File where predictions from glmnet will be saved.",
            metavar = "character"
        ),
        make_option(
            c("-l", "--lambda-output"),
            type    = "character",
            default = NULL,
            help    = "File where lambdas from glmnet will be saved.",
            metavar = "character"
        ),
        make_option(
            c("-b", "--beta-output"),
            type    = "character",
            default = NULL,
            help    = "File where betas (effect sizes) from glmnet will be saved.",
            metavar = "character"
        ),
        make_option(
            c("-B", "--BIM-file"),
            type    = "character",
            default = NULL,
            help    = "PLINK BIM file corresponding to current gene.",
            metavar = "character"
        ),
        make_option(
            c("-a", "--alpha"),
            type    = "double",
            default = 0.5,
            help    = "Alpha mixing parameter used for elastic net regression [default = %default].",
            metavar = "double"
        ),
        make_option(
            c("-n", "--num-folds"),
            type    = "integer",
            default = 10,
            help    = "Number of inner crossvalidation folds to use in the elastic net [default = %default].",
            metavar = "integer"
        ),
        make_option(
            c("-s", "--random-seed"),
            type    = "integer",
            default = 2018,
            help    = "Random seed used to produce reproducible crossvalidation fold splits [default = %default].",
            metavar = "integer"
        )
    )
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

    cat("Parsed options:\n\n")
    print(opt)

    # parse command-line arguments
    x.df.path = opt$genotype_dosage_file 
    y.df.path = opt$expression_file 
    gene.name = opt$gene_name 
    predout   = opt$prediction_output 
    lambdaout = opt$lambda_output 
    alpha     = opt$alpha 
    outpath   = opt$beta_output 
    bim.path  = opt$BIM_file 
    nfolds    = opt$num_folds 
    seed      = opt$random_seed

    # typecast numbers in case R read them as strings
    alpha  = as.numeric(alpha)
    nfolds = as.numeric(nfolds)

    # set random seed
    set.seed(seed)

    # ==========================================================================================
    # Directories & Variables
    # ==========================================================================================

    # which gene are we analyzing?
    cat(paste0("Building predictive model for gene ", gene.name, " with ", nfolds, "-fold crossvalidation.\n"))

    # load data
    # data frame x.df contains the genotype dosages for the current gene
    # data frame y.df contains the expression levels
    # nota bene: this script assumes that y.df is a data frame with column format
    # > "GENE,SAMPLE1,SAMPLE2,...,SAMPLEN"
    # containing inverse normal rank-norm'd RPKMs 
    # we parse data X and phenotype y from these data frames
    x.df = fread(x.df.path)
    y.df = fread(y.df.path, header = TRUE)

    # forcibly reorder rows of x.df
    setorder(x.df, IID)

    # use new column names to forcibly reorder y.df
    setorder(y.df, Gene)

    # now make a transpose
    ty = t(y.df[,-c(1)])
    colnames(ty) = y.df$Gene

    # add sample names to ty for reordering
    ty = ty[order(row.names(ty)),]

    # preallocate a vector to store predicted expression levels
    n = dim(ty)[1]
    y.pred = matrix(NA,n,1)

    # ensure that nfolds <= length(y) - 1
    # if not, then adjust and notify accordingly
    cat("Ensuring that nfolds <= length(y) - 1...\n")
    if ( nfolds <= n - 1) {
        cat("...looks good.\n")
    } else {
        cat("nfolds = ", nfolds, " but length(y) - 1 = ", n - 1, ", adjusting accordingly...\n")
        nfolds = n - 1
        cat("now nfolds = ", nfolds, "\n")
    }

    # prepare the gene output file that will contain the prediction and lambda for each held-out sample
    # put a header on it for bookkeeping
    # the pipeline will later discard the header when forming a unified prediction file, one line per prediction
    lambdaarray = c("Gene", "Held_out_sample", "Mean_MSE", "Lambda", "Prediction")
    write(lambdaarray, file = lambdaout, ncolumns = 5, sep = "\t", append = FALSE)

    # do same for weights (betas)
    betafile = c("Gene", "Held_out_sample", "SNP", "A1", "A2", "Beta")
    write(betafile, file = outpath, ncolumns = 6, append = FALSE, sep = "\t")

    for (held.out.sample in 1:n){

        # get NWDID of held-out sample
        held.out.sample.name = x.df$IID[held.out.sample]

        # we expect x.df to have a PLINK RAW format
        # thus, the FID/IIDs are the same,
        # and the dosages start at column 6
        # discard all columns 1-5 since we do not need them here
        tx        = data.frame(x.df[-held.out.sample, -c(1:6)])
        x.heldout = as.double(x.df[held.out.sample, -c(1:6)])
        names(x.heldout) = colnames(tx)

        y.heldout = ty[held.out.sample,] 
        ty.train  = ty[-held.out.sample,]

        # must "impute" missing dosages
        # use 2*MAF, since we don't need sophisticated dosage imputation
        # use only observed dosages; this means that means are not necessarily divided by full sample size!
        mafs.missing = 2 * apply(tx, 2, mean, na.rm = TRUE)
        for (i in seq_along(tx)){
            set(tx, j=i, value = as.numeric(tx[[i]]))
            set(tx, i = which( is.na(tx[[i]]) ), j=i, value = mafs.missing[i])
        }
        na.heldout = is.na(x.heldout)
        x.heldout[na.heldout] = mafs.missing[na.heldout]
        X = data.matrix(tx) # as opposed to as.matrix()...? see https://stackoverflow.com/questions/8458233/r-glmnet-as-matrix-error-message
        row.names(X) = x.df$IID[-held.out.sample]

        # pull SNPs with at least 1 minor allele
        # also forcibly reorder the rows (samples) of X
        # we did this for the expression levels ty too
        minorsnps = subset(colMeans(X), colMeans(X, na.rm=TRUE)>0)
        minorsnps = names(minorsnps)
        X = X[,minorsnps]

        # ensure that X and y have the same rows
        matching.samples = row.names(ty.train) %in% row.names(X)
        ty.train = subset(ty.train, matching.samples)
        X = X[matching.samples,]

        # mean-center the genotype variables
        X = scale(X, center = TRUE, scale = FALSE)

        # initialize data.frame of best betas
        bestbetas = data.frame()

        # predicted expression for held-out sample is missing until calculated otherwise
        y.heldout.expr.pred = NA

        # skip genes with < 2 cis-SNPs
        if(is.null(dim(X)) | dim(X)[2] == 0){
            bestbetas = data.frame() ###effectively skips genes with <2 cis-SNPs
        } else {
            # extract phenotype from expression matrix
            # the grep command pulls the numerical position of the gene in the matrix
            # save that integer and use for subsetting
            y.idx = grep(gene.name, y.df$Gene) 

            if (length(y.idx) != 0){

                # this makes a numeric matrix for the phenotype
                # all missing quantities are first set to 0
                y = as.vector(ty.train[,gene.name])
                y[is.na(y)] = 0

                if (length(y) > 0){

                    # want to prevent unusual error caused by random generation of lambda values:
                    # https://github.com/lmweber/glmnet-error-example/blob/master/glmnet_error_example.R
                    # therefore, find lambda.max per glmnet default and set lambdas manually
                    sx   = as.matrix(scale(X, center = TRUE, scale = TRUE))
                    sy   = as.vector(scale(y, center = TRUE, scale = TRUE))
                    lam.max = max(abs(colSums(sx*sy))) / n
                    nlambda = 100
                    lambda  = exp(seq(log(0.01), log(lam.max), length.out=nlambda))

                    # fix cross-validation folds for reproducibility
                    # note: this makes "nfolds" argument to cv.glmnet redundant
                    flds    = createFolds(y, k = nfolds, list = TRUE, returnTrain = FALSE)
                    foldids = rep(1, length(y))
                    for (i in 1:nfolds) {
                        foldids[flds[[i]]] = i
                    }

                    # cross-validate!
                    cat(paste0("crossvalidating sample ", held.out.sample.name, "..."))
                    #my.glmnet = cv.glmnet(x=X, y=y, nfolds=nfolds, family = "gaussian", alpha=alpha, keep=TRUE, grouped = FALSE, dfmax=n, pmax=n, lambda.min.ratio = 0.01, nlambda=nlambda)
                    my.glmnet = cv.glmnet(x=X, y=y, nfolds=nfolds, family = "gaussian", alpha=alpha, keep=TRUE, grouped = FALSE, dfmax=n, pmax=n, lambda.min.ratio = 0.01, nlambda=nlambda, foldid=foldids)
                    cat(paste0("done at timestamp ", Sys.time(), "\n"))

                    # pull info to find best lambda
                    fit.df = data.frame("cvm" = my.glmnet$cvm, "lambda" = my.glmnet$lambda, "nrow" = 1:length(my.glmnet$cvm)) 

                    # parse results
                    best.lam    = fit.df[which.min(fit.df[,1]),]                         # use which.min or which.max depending on cv measure (MSE min, AUC max, ...)
                    cvm.best    = best.lam[,1]                                           # best crossvalidated mean squared error
                    lambda.best = best.lam[,2]                                           # corresponding λ for best MSE
                    nrow.best   = best.lam[,3]                                           # position of best λ in cv.glmnet output
                    ret         = as.data.frame(my.glmnet$glmnet.fit$beta[,nrow.best])   # get βs from best λ 
                    ret[ret == 0.0] = NA
                    bestbetas   = as.vector(ret[which(!is.na(ret)),])                    # vector of nonzero βs 
                    names(bestbetas) = rownames(ret)[which(!is.na(ret))]
                    pred.mat    = as.matrix(my.glmnet$fit.preval[,nrow.best])            # pull out predictions at best λ 

                    # make prediction on held-out sample
                    # remember to only use SNPs that met minor allele threshold
                    # can use s = lambda.min to use best predictive λ directly
                    y.heldout.expr.pred = predict(my.glmnet, newx=t(as.matrix(x.heldout[minorsnps])), s = my.glmnet$lambda.min)

                    # save information about each internal CV fold
                    # order is gene, NWDID of held out sample, cvm of internal LOOCV, corresponding lambda 
                    lambdaarray = c(gene.name, held.out.sample.name, my.glmnet$cvm[which.min(my.glmnet$cvm)], my.glmnet$lambda.min, y.heldout.expr.pred)
                    write(lambdaarray, file = lambdaout, ncolumns = 5, sep = "\t", append = TRUE)

                    # record results when glmnet finds an eQTL and that eQTL is in our genotype set
                    # otherwise record missing values
                    betafile = c(gene.name, held.out.sample.name, NA, NA, NA, NA)
                    if(length(bestbetas) > 0 & !is.null(dim(pred.mat))){
                        # output best shrunken betas for PrediXcan
                        # output format: "gene","SNP","refAllele","effectAllele","beta"
                        # note that this makes explicit reference to rsid format from PLINK BIM from SAGE merged LATplus array
                        bestbetalist = names(bestbetas) # next lines format entries of this list to match BIM
                        bestbetalist = gsub("X", "", bestbetalist)
                        bestbetalist = gsub(".", ":", bestbetalist, fixed = TRUE)
                        bestbetalist = gsub("_.*", "", bestbetalist, perl = TRUE)
                        bimfile      = fread(bim.path)
                        bimfile$V2   = gsub("-", ":", bimfile$V2) # needed for SAGE LAT-LAT+ merged genotypes
                        bestbetainfo = bimfile[bimfile$V2 %in% bestbetalist,]
                        betatable    = as.matrix(cbind(bestbetainfo, bestbetas))
                        betafile     = cbind(gene.name, held.out.sample.name, betatable[,2], betatable[,5], betatable[,6], betatable[,7])
                    }

                    # save betas to file
                    # t() necessary for correct output from write() function
                    write(t(betafile), file = outpath, ncolumns = 6, append = TRUE, sep = "\t")

                } else{
                    cat(paste("no expression data for", gene.name, "\nMarking missing results for ", x.df$IID[held.out.sample], "\n"))
                    lambdaarray = c(gene.name, x.df$IID[held.out.sample], NA, -1.0, NA) 
                    write(lambdaarray, file = lambdaout, ncolumns = 5, sep = "\t", append = TRUE)
                }
            } else{
                cat(paste("no expression data for", gene.name, "\nMarking missing results for ", x.df$IID[held.out.sample], "\n"))
                lambdaarray = c(gene.name, x.df$IID[held.out.sample], NA, -1.0, NA) 
                write(lambdaarray, file = lambdaout, ncolumns = 5, sep = "\t", append = TRUE)
            }

        # save predicted expression level
        y.pred[held.out.sample] = y.heldout.expr.pred

        }

    }

    # write y.pred to file
    predictionarray = c(gene.name, t(y.pred))
    write(predictionarray, file = predout, ncolumns = n + 1, sep = "\t", append = FALSE)

    return()
}

# ==========================================================================================
# Executable script code
# ==========================================================================================

# timestamp
cat(paste("Begin estimating prediction weights at date and time ", Sys.time(), "\n", sep = "") )

# run function
compute.new.weights()

# any warnings?
# note that glmnet may complain whenever its Fortran routine executes early on the lambda path
# that is ok
# if it complains about missing Fortran routines, then take them seriously
# the latter warnings indicate a missing glmnet installation, a mismatched Rscript/glmnet version,
# or an unexpected version of Rscript (e.g. 3.5 vs. 3.4)
cat("any warnings?\n")
warnings()

cat(paste("Done estimating prediction weights at date and time ", Sys.time(), "\n", sep = "") )
