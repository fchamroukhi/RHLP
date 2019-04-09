rm(list = ls())
source("R/FData.R")
source("R/ModelRHLP.R")
source("R/enums.R")
source("R/EM.R")

fileName = "data/simulated_time_series.mat"
library(R.matlab)
data <- readMat(fileName)
#construction des matrices de regression
X <- data$x
Y <- data$y
fData <- FData$new()
fData$setData(X, Y)


K <- 5; # nomber of regimes (mixture components)
p <- 3; # dimension of beta' (order of the polynomial regressors)
q <- 1; # dimension of w (ordre of the logistic regression: to be set to 1 for segmentation)
variance_type <- variance_types$hetereskedastic

n_tries=1
max_iter=1500
threshold <- 1e-6
verbose <- TRUE
verbose_IRLS <- FALSE

####
# EM Algorithm
####
solution <- EM(K, p, q, variance_type, fData, n_tries, max_iter, threshold, verbose, verbose_IRLS) # instance of class StatRHLP

solution$paramSolution
solution$statSolution

