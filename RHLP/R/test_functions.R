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
