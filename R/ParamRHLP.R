#' A Reference Class which contains parameters of a RHLP model.
#'
#' ParamRHLP contains all the parameters of a [RHLP][ModelRHLP] model.
#'
#' @usage NULL
#' @field W Parameters of the logistic process.
#' \eqn{W = w_{1},\dots,w_{K-1}}{W = (w1,\dots,wK-1)} is a matrix of dimension
#' \eqn{(q + 1, K - 1)}, with \emph{q} the order of the logistic
#' regression.
#' @field beta Parameters of the polynomial regressions.
#' \eqn{\beta = (\beta_{1},\dots,\beta_{K})}{\beta = (\beta1,\dots,\betaK)} is
#' a matrix of dimension \eqn{(p + 1, K)}, with \emph{p} the
#' order of the polynomial regression.
#' @field sigma The variances for the \emph{K} regimes. If [RHLP][ModelRHLP]
#' model is homoskedastic (\emph{variance_type} = 1) then sigma is a matrix of
#' size \eqn{(1, 1)}, else if [RHLP][ModelRHLP] model is heteroskedastic then
#' sigma is a matrix of size \eqn{(K, 1)}.
#' @seealso [ModelRHLP]
#' @export
ParamRHLP <- setRefClass(
  "ParamRHLP",
  fields = list(
    W = "matrix",
    beta = "matrix",
    sigma = "matrix"
  ),
  methods = list(
    initParam = function(modelRHLP, phi, try_algo = 1) {

      if (try_algo == 1) { # Uniform segmentation into K contiguous segments, and then a regression

        # Initialization of W
        W <<- zeros(modelRHLP$q + 1, modelRHLP$K - 1)

        zi <- round(modelRHLP$m / modelRHLP$K) - 1

        beta <<- matrix(NA, modelRHLP$p + 1, modelRHLP$K)

        for (k in 1:modelRHLP$K) {
          i <- (k - 1) * zi + 1
          j <- k * zi
          yij <- modelRHLP$Y[i:j]

          Phi_ij <- phi$XBeta[i:j, ]

          bk <-  solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% yij
          beta[, k] <<- bk

          if (modelRHLP$variance_type == variance_types$homoskedastic) {
            sigma <<- matrix(1)
          }
          else{
            sigma[k] <<- var(yij)
          }
        }
      } else{# Random segmentation into K contiguous segments, and then a regression

        # Initialization of W
        W <<- rand(modelRHLP$q + 1, modelRHLP$K - 1)

        Lmin <- round(modelRHLP$m / (modelRHLP$K + 1)) # Minimum number of points in a segment
        tk_init <- zeros(modelRHLP$K, 1)
        if (modelRHLP$K == 1) {
          tk_init[2] = modelRHLP$m
        } else {
          K_1 <- modelRHLP$K
          for (k in 2:modelRHLP$K) {
            K_1 <- (K_1 - 1)
            temp <-
              (tk_init[k - 1] + Lmin):(modelRHLP$m - (K_1 * Lmin))

            ind <- sample(length(temp))

            tk_init[k] <- temp[ind[1]]
          }
          tk_init[modelRHLP$K + 1] <- modelRHLP$m
        }

        beta <<- matrix(NA, modelRHLP$p + 1, modelRHLP$K)
        for (k in 1:modelRHLP$K) {
          i <- tk_init[k] + 1
          j <- tk_init[k + 1]
          yij <- modelRHLP$Y[i:j]
          Phi_ij <- phi$XBeta[i:j, ]
          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% yij
          beta[, k] <<- bk

          if (modelRHLP$variance_type == variance_types$homoskedastic) {
            sigma <<- var(modelRHLP$Y)
          }
          else{
            sigma[k] <<- 1
          }
        }
      }
    },

    MStep = function(modelRHLP, statRHLP, phi, verbose_IRLS) {

      # Maximization w.r.t betak and sigmak (the variances)
      if (modelRHLP$variance_type == variance_types$homoskedastic) {
        s = 0
      }
      for (k in 1:modelRHLP$K) {
        weights <- statRHLP$tik[, k] # Post prob of each component k (dimension nx1)
        nk <- sum(weights) # Expected cardinal number of class k

        Xk <- phi$XBeta * (sqrt(weights) %*% ones(1, modelRHLP$p + 1)) # [m*(p+1)]
        yk <- modelRHLP$Y * (sqrt(weights)) # Dimension: (nx1).*(nx1) = (nx1)

        M <- t(Xk) %*% Xk
        epps <- 1e-9
        M <- M + epps * diag(modelRHLP$p + 1)

        beta[, k] <<- solve(M) %*% t(Xk) %*% yk # Maximization w.r.t betak
        z <- sqrt(weights) * (modelRHLP$Y - phi$XBeta %*% beta[, k])

        # Maximisation w.r.t sigmak (the variances)
        sk <- t(z) %*% z

        if (modelRHLP$variance_type == variance_types$homoskedastic) {
          s <- s + sk
          sigma <<- s / modelRHLP$m
        } else {
          sigma[k] <<- sk / nk
        }
      }

      # Maximization w.r.t W
      #  IRLS : Iteratively Reweighted Least Squares (for IRLS, see the IJCNN 2009 paper)
      res_irls <- IRLS(phi$Xw, statRHLP$tik, ones(nrow(statRHLP$tik), 1), W, verbose_IRLS)

      W <<- res_irls$W
      piik <- res_irls$piik
      reg_irls <- res_irls$reg_irls
    }
  )
)

ParamRHLP <- function(modelRHLP) {
  W <- matrix(0, modelRHLP$p + 1, modelRHLP$K - 1)
  beta <- matrix(NA, modelRHLP$p + 1, modelRHLP$K)
  if (modelRHLP$variance_type == variance_types$homoskedastic) {
    sigma <- matrix(NA)
  }
  else{
    sigma <- matrix(NA, modelRHLP$K)
  }
  new("ParamRHLP", W = W, beta = beta, sigma = sigma)
}
