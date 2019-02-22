# in constructing PrediXcan weights, the PrediXcan team did not standardize the genotypes before running
# see here:
# https://github.com/hakyimlab/PrediXcan/blob/e69766838a71c139bf8bb7d745b8938365e2d029/Paper-Scripts/Heather/DGN-calc-weights/01_imputedDGN-WB_CV_elasticNet.r#L116
# we mean-center the genotypes before running
# this nuance should not change results
# here is a minimal working example

# load libraries
library(glmnet)
library(caret)

# get a toy dataset
data(iris)
iris = iris
dat = iris[iris$Species %in% c("setosa","versicolor"),]
X = as.matrix(dat[,1:4])
Y = as.factor(as.character(dat$Species))
my.seed = 284
lambdas = rev( seq(0, 1, length = 100) )

# make reproducible folds
set.seed(my.seed)
flds <- createFolds(Y, k = 5, list = TRUE, returnTrain = FALSE)
foldids = rep(1,length(Y))
foldids[flds$Fold2] = 2
foldids[flds$Fold3] = 3
foldids[flds$Fold4] = 4
foldids[flds$Fold5] = 5

# don't use glmnet standardize
set.seed(my.seed)
model2 = cv.glmnet(x = X, y = Y, family = "binomial", standardize = FALSE, alpha = 0.5, lambda = lambdas, nfolds = 5, foldid = foldids)

# with glmnet standardize
set.seed(my.seed)
model4 = cv.glmnet(x = scale(X, center = T, scale = F), y = Y, family = "binomial", standardize = TRUE, alpha = 0.5, lambda = lambdas, nfolds = 5, foldid = foldids)

# verify that coefficients are exactly the same
sum(model4$glmnet.fit$beta - model2$glmnet.fit$beta) ## should be 0
