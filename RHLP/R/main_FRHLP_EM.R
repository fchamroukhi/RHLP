rm(list = ls())
source("R/dataset.R")
source("R/MixModel.R")
source("R/ModelOptions.R")
source("R/enums.R")
source("R/ModelLearner.R")

fileName = "data/simulated_time_series.mat"
mixData <- MyData$new()
mixData$setDataFromMat(fileName)
K <- 5; # nomber of regimes (mixture components)
p <- 3; # dimension of beta' (order of the polynomial regressors)
q <- 1; # dimension of w (ordre of the logistic regression: to be set to 1 for segmentation)
mixModel <- MixModel(mixData,K,p,q)

n_tries=1
max_iter=1500
threshold <- 1e-6
verbose <- TRUE
verbose_IRLS <- FALSE
modelOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, variance_types$hetereskedastic)

####
# EM Algorithm
####

solution <- EM(mixModel, modelOptions)
#mixParamSolution <- solution[[1]]
#mixStatsSolution <- solution[[2]]

par(mfrow=c(2, 1))
plot(mixData$y, type = "l", ylab = "y", xlab = "")
title(main = "Time series, RHLP regimes and process probabilities")
colors = rainbow(K)
for (k in 1:K) {
  index <- (solution[[2]]$klas == k)
  polynomials <- solution[[2]]$polynomials[(solution[[2]]$klas == k), k]
  lines(solution[[2]]$polynomials[, k], lty = "dotted", lwd = 2, col = colors[k])
  lines(seq(1:mixData$m)[index], polynomials, lwd = 2, col = colors[k])
}

plot(solution[[2]]$h_ig[, 1], type = "l", lwd = 2, col = colors[1], xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)))
for (k in 2:K) {
  lines(solution[[2]]$h_ig[, k], type = "l", lwd = 2, col = colors[k])
}

# plot
par(mfrow=c(2, 1))
plot(mixData$y, type = "l", ylab = "y", xlab = "", col = "grey")
title(main = "Time series, estimated RHLP model, and segmentation")

tk = which(diff(solution[[2]]$klas) != 0)
for (i in 1:length(tk)) {
  abline(v = tk[i], lty = "dotted", lwd = 2, col = "black")
}
for (k in 1:K) {
  index <- (solution[[2]]$klas == k)
  polynomials <- solution[[2]]$polynomials[(solution[[2]]$klas == k), k]
  lines(seq(1:mixData$m)[index], polynomials, lwd = 2, col = "red")
}
plot(solution[[2]]$klas, type = "l", lwd = 2, col = "black", xlab = "", ylab = "Estimated class labels")
