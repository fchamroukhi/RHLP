source("R/dataset.R")
source("R/enums.R")
source("R/ModelOptions.R")

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
