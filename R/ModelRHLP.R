#' A Reference Class which represents a fitted RHLP model.
#'
#' ModelRHLP represents a [RHLP][ModelRHLP] model for which parameters have
#' been estimated.
#'
#' @usage NULL
#' @field paramRHLP A [ParamRHLP][ParamRHLP] object. It contains the estimated values of the parameters.
#' @field statRHLP A [StatRHLP][StatRHLP] object. It contains all the statistics associated to the RHLP model.
#' @seealso [ParamRHLP], [StatRHLP]
#' @export
ModelRHLP <- setRefClass(
  "ModelRHLP",
  fields = list(
    paramRHLP = "ParamRHLP",
    statRHLP = "StatRHLP"
  ),
  methods = list(

    plot = function() {

      oldpar <- par()[c("mfrow", "mai", "mgp")]
      on.exit(par(oldpar), add = TRUE)

      yaxislim <- c(mean(paramRHLP$fData$Y) - 2 * sd(paramRHLP$fData$Y), mean(paramRHLP$fData$Y) + 2 * sd(paramRHLP$fData$Y))

      # Data, regressors, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(paramRHLP$fData$X, paramRHLP$fData$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      title(main = "Time series, RHLP regimes and process probabilities")
      colorsvec = rainbow(paramRHLP$K)
      for (k in 1:paramRHLP$K) {
        index <- (statRHLP$klas == k)
        polynomials <- statRHLP$polynomials[index, k]
        lines(paramRHLP$fData$X, statRHLP$polynomials[, k], col = colorsvec[k], lty = "dotted", lwd = 1.5)
        lines(paramRHLP$fData$X[index], col = colorsvec[k], polynomials, lwd = 1.5)
      }

      # Probablities of the hidden process (segmentation)
      plot.default(paramRHLP$fData$X, statRHLP$pi_ik[, 1], type = "l", xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)), col = colorsvec[1], lwd = 1.5, ylim = c(0, 1))
      if (paramRHLP$K > 1) {
        for (k in 2:paramRHLP$K) {
          lines(paramRHLP$fData$X, statRHLP$pi_ik[, k], col = colorsvec[k], lwd = 1.5, ylim = c(0, 1))
        }
      }

      # Data, regression model, and segmentation
      par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
      plot.default(paramRHLP$fData$X, paramRHLP$fData$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
      lines(paramRHLP$fData$X, statRHLP$Ex, col = "red", lwd = 1.5)
      title(main = "Time series, estimated RHLP model, and segmentation")

      # Transition time points
      tk <- which(diff(statRHLP$klas) != 0)
      for (i in 1:length(tk)) {
        abline(v = tk[i], col = "red", lty = "dotted", lwd = 1.5)
      }

      # Probablities of the hidden process (segmentation)
      plot.default(paramRHLP$fData$X, statRHLP$klas, type = "l", xlab = "x", ylab = "Estimated class labels", col = "red", lwd = 1.5, yaxt = "n")
      axis(side = 2, at=1:paramRHLP$K)
    }
  )
)
