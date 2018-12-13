rm(list = ls())
source("R/dataset.R")
source("R/MixModel.R")
source("R/ModelOptions.R")
source("R/enums.R")
source("R/ModelLearner.R")

library(R.matlab)
fileName = "data/simulated_time_series.mat"
mixData <- MyData$new()
mixData$setDataFromMat(fileName)
K <- 5; # nomber of regimes (mixture components)
p <- 3; # dimension of beta' (order of the polynomial regressors)
q <- 1; # dimension of w (ordre of the logistic regression: to be set to 1 for segmentation)
mixModel <- MixModel(mixData,K,p,q)

n_tries=1
max_iter=1000
threshold <- 1e-5
verbose <- TRUE
verbose_IRLS <- FALSE
modelOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, variance_types$homoskedastic)

####
# EM Algorithm
####

solution <- EM(mixModel, modelOptions)
#mixParamSolution <- solution[[1]]
#mixStatsSolution <- solution[[2]]

