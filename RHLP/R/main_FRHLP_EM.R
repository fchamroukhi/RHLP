rm(list = ls())
source("R/FData.R")
source("R/ModelRHLP.R")
source("R/enums.R")
source("R/ModelLearner.R")


# Building matrices for regression
load("data/simulatedTimeSeries.RData")
fData <- FData$new()
fData$setData(X, Y)


K <- 5 # number of regimes (mixture components)
p <- 3 # dimension of beta (order of the polynomial regressors)
q <- 1 # dimension of w (order of the logistic regression: to be set to 1 for segmentation)
variance_type <- variance_types$hetereskedastic

modelRHLP <- ModelRHLP(fData, K, p, q, variance_type)

n_tries <- 1
max_iter = 1500
threshold <- 1e-6
verbose <- TRUE
verbose_IRLS <- FALSE

####
# EM Algorithm
####
solution <- EM(modelRHLP, n_tries, max_iter, threshold, verbose, verbose_IRLS)

solution$plot()
