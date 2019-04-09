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


par(mfrow=c(2, 1))
plot(fData$Y, type = "l", ylab = "y", xlab = "")
title(main = "Time series, RHLP regimes and process probabilities")
colors = rainbow(solution$paramSolution$K)
for (k in 1:solution$paramSolution$K) {
  index <- (solution$statSolution$klas == k)
  polynomials <- solution$statSolution$polynomials[(solution$statSolution$klas == k), k]
  lines(solution$statSolution$polynomials[, k], lty = "dotted", lwd = 2, col = colors[k])
  lines(seq(1:fData$m)[index], polynomials, lwd = 2, col = colors[k])
}

plot(solution$statSolution$piik[, 1], type = "l", lwd = 2, col = colors[1], xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)))
for (k in 2:solution$paramSolution$K) {
  lines(solution$statSolution$piik[, k], type = "l", lwd = 2, col = colors[k])
}

# plot
par(mfrow=c(2, 1))
plot(fData$Y, type = "l", ylab = "y", xlab = "")
title(main = "Time series, estimated RHLP model, and segmentation")

tk = which(diff(solution$statSolution$klas) != 0)
for (i in 1:length(tk)) {
  abline(v = tk[i], lty = "dotted", lwd = 2, col = "red")
}
plot(solution$statSolution$klas, type = "l", lwd = 2, col = "red", xlab = "", ylab = "Estimated class labels")

