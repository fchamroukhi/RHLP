source("R/dataset.R")
testLoadDataFromMat <- function(){
  library(R.matlab)
  fileName = "data/simulated_time_series.mat"
  data <- MyData$new()
  data$setDataFromMat(fileName)
  #x <- data$x
  #y <- data$y
  print(data$m)
}
testLoadDataFromMat()

n_tries=1
max_iter=1000
threshold <- 1e-5
verbose <- TRUE
verbose_IRLS <- TRUE
mixOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, variance_types$common)
