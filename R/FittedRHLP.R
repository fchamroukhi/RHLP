FittedRHLP <- setRefClass(
  "FittedRHLP",
  fields = list(
    modelRHLP = "ModelRHLP",
    paramRHLP = "ParamRHLP",
    statRHLP = "StatRHLP"
  ),
  methods = list(
    plot = function() {

      oldpar <- par()[c("mfrow", "mai", "mgp")]
      on.exit(par(oldpar), add = TRUE)

      yaxislim <- c(mean(modelRHLP$Y) - 2 * sd(modelRHLP$Y), mean(modelRHLP$Y) + 2 * sd(modelRHLP$Y))

      # Data, regressors, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(modelRHLP$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      title(main = "Time series, RHLP regimes and process probabilities")
      colorsvec = rainbow(modelRHLP$K)
      for (k in 1:modelRHLP$K) {
        index <- (statRHLP$klas == k)
        polynomials <- statRHLP$polynomials[index, k]
        lines(statRHLP$polynomials[, k], col = colorsvec[k], lty = "dotted", lwd = 1.5)
        lines(seq(1:modelRHLP$m)[index], col = colorsvec[k], polynomials, lwd = 1.5)
      }

      # Probablities of the hidden process (segmentation)
      plot.default(statRHLP$piik[, 1], type = "l", xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)), col = colorsvec[1], lwd = 1.5)
      if (modelRHLP$K > 1) {
        for (k in 2:modelRHLP$K) {
          lines(statRHLP$piik[, k], col = colorsvec[k], lwd = 1.5)
        }
      }

      # Data, regression model, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(modelRHLP$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      lines(statRHLP$Ex, col = "red", lwd = 1.5)
      title(main = "Time series, estimated RHLP model, and segmentation")

      # Transition time points
      tk <- which(diff(statRHLP$klas) != 0)
      for (i in 1:length(tk)) {
        abline(v = tk[i], col = "red", lty = "dotted", lwd = 1.5)
      }

      # Probablities of the hidden process (segmentation)
      plot.default(statRHLP$klas, type = "l", xlab = "x", ylab = "Estimated class labels", col = "red", lwd = 1.5)
    }
  )
)

FittedRHLP <- function(modelRHLP, paramRHLP, statRHLP) {
  new("FittedRHLP", modelRHLP = modelRHLP, paramRHLP = paramRHLP, statRHLP = statRHLP)
}
