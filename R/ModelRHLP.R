#' A Reference Class which represents a fitted RHLP model.
#'
#' ModelRHLP represents an estimated RHLP model.
#'
#' @field paramRHLP A [ParamRHLP][ParamRHLP] object. It contains the estimated
#'   values of the parameters.
#' @field statRHLP A [StatRHLP][StatRHLP] object. It contains all the
#'   statistics associated to the RHLP model.
#' @seealso [ParamRHLP], [StatRHLP]
#' @export
ModelRHLP <- setRefClass(
  "ModelRHLP",
  fields = list(
    paramRHLP = "ParamRHLP",
    statRHLP = "StatRHLP"
  ),
  methods = list(

    plot = function(what = c("regressors", "estimatedsignal")) {
      "Plot method.
      \\describe{
        \\item{\\code{what}}{The type of graph requested:
          \\itemize{
            \\item \\code{\"regressors\" = } Polynomial regression components.
            \\item \\code{\"estimatedsignal\" = } Estimated signal (Sum of the
            polynomial components weighted by the logistic probabilities.
          }
        }
      }
      By default, all the above graphs are produced."

      what <- match.arg(what, several.ok = TRUE)

      oldpar <- par()[c("mfrow", "mai", "mgp")]
      on.exit(par(oldpar), add = TRUE)

      yaxislim <- c(mean(paramRHLP$Y) - 2 * sd(paramRHLP$Y), mean(paramRHLP$Y) + 2 * sd(paramRHLP$Y))

      if (any(what == "regressors")) {
        # Data, regressors, and segmentation
        par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
        plot.default(paramRHLP$X, paramRHLP$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
        title(main = "Time series, RHLP regimes and process probabilities")
        colorsvec = rainbow(paramRHLP$K)
        for (k in 1:paramRHLP$K) {
          index <- (statRHLP$klas == k)
          polynomials <- statRHLP$polynomials[index, k]
          lines(paramRHLP$X, statRHLP$polynomials[, k], col = colorsvec[k], lty = "dotted", lwd = 1.5)
          lines(paramRHLP$X[index], col = colorsvec[k], polynomials, lwd = 1.5)
        }

        # Probablities of the hidden process (segmentation)
        plot.default(paramRHLP$X, statRHLP$pi_ik[, 1], type = "l", xlab = "x", ylab = expression('Probability ' ~ pi [k] (t, w)), col = colorsvec[1], lwd = 1.5, ylim = c(0, 1))
        if (paramRHLP$K > 1) {
          for (k in 2:paramRHLP$K) {
            lines(paramRHLP$X, statRHLP$pi_ik[, k], col = colorsvec[k], lwd = 1.5, ylim = c(0, 1))
          }
        }
      }

      if (any(what == "estimatedsignal")) {
        # Data, regression model, and segmentation
        par(mfrow = c(2, 1), mai = c(0.6, 1, 0.5, 0.5), mgp = c(2, 1, 0))
        plot.default(paramRHLP$X, paramRHLP$Y, type = "l", ylim = yaxislim, xlab = "x", ylab = "y")
        lines(paramRHLP$X, statRHLP$Ex, col = "red", lwd = 1.5)
        title(main = "Time series, estimated RHLP model, and segmentation")

        # Transition time points
        tk <- which(diff(statRHLP$klas) != 0)
        for (i in 1:length(tk)) {
          abline(v = paramRHLP$X[tk[i]], col = "red", lty = "dotted", lwd = 1.5)
        }

        # Probablities of the hidden process (segmentation)
        plot.default(paramRHLP$X, statRHLP$klas, type = "l", xlab = "x", ylab = "Estimated class labels", col = "red", lwd = 1.5, yaxt = "n")
        axis(side = 2, at = 1:paramRHLP$K)
      }
    },

    summary = function() {
      "Summary method."
      digits = getOption("digits")

      title <- paste("Fitted RHLP model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(paste0("RHLP model with K = ", paramRHLP$K, ifelse(paramRHLP$K > 1, " components", " component"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = statRHLP$loglik, "nu" = paramRHLP$nu, "AIC" = statRHLP$AIC,
                        "BIC" = statRHLP$BIC, "ICL" = statRHLP$ICL, row.names = "", check.names = FALSE)
      print(tab, digits = digits)

      cat("\nClustering table (Number of observations in each regimes):\n")
      print(table(statRHLP$klas))

      cat("\nRegression coefficients:\n\n")
      if (paramRHLP$p > 0) {
        row.names = c("1", sapply(1:paramRHLP$p, function(x) paste0("X^", x)))
      } else {
        row.names = "1"
      }

      betas <- data.frame(paramRHLP$beta, row.names = row.names)
      colnames(betas) <- sapply(1:paramRHLP$K, function(x) paste0("Beta(K = ", x, ")"))
      print(betas, digits = digits)

      cat(paste0(ifelse(paramRHLP$variance_type == "homoskedastic", "\n\n",
                        "\nVariances:\n\n")))
      sigma2 = data.frame(t(paramRHLP$sigma2), row.names = NULL)
      if (paramRHLP$variance_type == "homoskedastic") {
        colnames(sigma2) = "Sigma2"
        print(sigma2, digits = digits, row.names = FALSE)
      } else {
        colnames(sigma2) = sapply(1:paramRHLP$K, function(x) paste0("Sigma2(K = ", x, ")"))
        print(sigma2, digits = digits, row.names = FALSE)
      }

    }
  )
)
