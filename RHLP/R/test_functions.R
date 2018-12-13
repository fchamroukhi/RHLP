rm(list = ls())
source("R/dataset.R")
source("R/enums.R")
source("R/ModelOptions.R")
source("R/MixModel.R")
source("R/Phi.R")

testLoadDataFromMat <- function(){
  library(R.matlab)
  fileName = "data/simulated_time_series.mat"
  data <- MyData$new()
  data$setDataFromMat(fileName)
  #print(data$m)
}
testLoadDataFromMat()

testModelOptions <- function(){
  n_tries=1
  max_iter=1000
  threshold <- 1e-5
  verbose <- TRUE
  verbose_IRLS <- TRUE
  mixOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, variance_types$homoskedastic)
}
testModelOptions()

testModel <- function(){
  library(R.matlab)
  fileName = "data/simulated_time_series.mat"
  data <- MyData$new()
  data$setDataFromMat(fileName)

  K <- 3; #nombre de regimes
  p <- 1; #dimension de beta (ordre de reg polynomiale)
  q <- 1; #dimension de w (ordre de reg logistique)
  mix <- MixModel(data,K,p,q)
}
testModel()

testPhi <- function(){
  library(R.matlab)
  fileName = "data/simulated_time_series.mat"
  data <- MyData$new()
  data$setDataFromMat(fileName)

  K <- 3; #nombre de regimes
  p <- 1; #dimension de beta (ordre de reg polynomiale)
  q <- 1; #dimension de w (ordre de reg logistique)
  model <- MixModel(data,K,p,q)

  phi <- Phi$new()
  phi$setPhi1(model$x,model$p,model$q)
  return(phi)
}
phi <- testPhi()


testParameterInitialization <- function(){
  library(R.matlab)
  fileName = "data/simulated_time_series.mat"
  data <- MyData$new()
  data$setDataFromMat(fileName)

  K <- 3; #nombre de regimes
  p <- 1; #dimension de beta (ordre de reg polynomiale)
  q <- 1; #dimension de w (ordre de reg logistique)
  mix <- MixModel(data,K,p,q)

  n_tries=1
  max_iter=1000
  threshold <- 1e-5
  verbose <- TRUE
  verbose_IRLS <- TRUE
  mixOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, variance_types$homoskedastic)

  phi <- Phi$new()
  phi$setPhi1(mix$x,mix$p,mix$q)

  param <- MixParam(mix, mixOptions)
  param$initParam(mix, phi, mixOptions, try_algo = 1)

  return(param)
}
param <- testParameterInitialization()
