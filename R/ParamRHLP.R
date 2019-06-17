#' A Reference Class which contains parameters of a RHLP model.
#'
#' ParamRHLP contains all the parameters of a RHLP model.
#'
#' @field fData [FData][FData] object representing the sample.
#' @field K The number of regimes (mixture components).
#' @field p The order of the polynomial regression.
#' @field q The dimension of the logistic regression. For the purpose of
#' segmentation, it must be set to 1.
#' @field variance_type Numeric indicating if the model is homoskedastic
#' (`variance_type` = 1) or heteroskedastic (`variance_type` = 2).
#' @field W Parameters of the logistic process.
#' \eqn{W = w_{1},\dots,w_{K-1}}{W = (w1,\dots,wK-1)} is a matrix of dimension
#' \eqn{(q + 1, K - 1)}, with \emph{q} the order of the logistic regression.
#' @field beta Parameters of the polynomial regressions.
#' \eqn{\beta = (\beta_{1},\dots,\beta_{K})}{\beta = (\beta1,\dots,\betaK)} is
#' a matrix of dimension \eqn{(p + 1, K)}, with \emph{p} the order of the
#' polynomial regression.
#' @field sigma2 The variances for the \emph{K} regimes. If RHLP model is
#' homoskedastic (\emph{variance_type} = 1) then sigma2 is a matrix of size
#' \eqn{(1, 1)}, else if RHLP model is heteroskedastic then sigma2 is a matrix
#' of size \eqn{(K, 1)}.
#' @seealso [FData]
#' @export
ParamRHLP <- setRefClass(
  "ParamRHLP",
  fields = list(
    fData = "FData",
    phi = "list",

    K = "numeric",
    p = "numeric",
    q = "numeric",
    variance_type = "numeric",
    nu = "numeric",

    W = "matrix",
    beta = "matrix",
    sigma2 = "matrix"
  ),
  methods = list(

    initialize = function(fData = FData(numeric(1), matrix(1)), K = 1, p = 2, q = 1, variance_type = 1) {

      fData <<- fData

      phi <<- designmatrix(x = fData$X, p = p, q = q)

      K <<- K
      p <<- p
      q <<- q
      variance_type <<- variance_type

      if (variance_type == variance_types$homoskedastic) {
        nu <<- (p + q + 3) * K - (q + 1) - (K - 1)
      } else{
        nu <<- (p + q + 3) * K - (q + 1)
      }

      W <<- matrix(0, p + 1, K - 1)
      beta <<- matrix(NA, p + 1, K)

      if (variance_type == variance_types$homoskedastic) {
        sigma2 <<- matrix(NA)
      }
      else{
        sigma2 <<- matrix(NA, K)
      }

    },

    initParam = function(try_algo = 1) {
      "Method to initialize parameters \\code{W}, \\code{beta} and
      \\code{sigma2}.

      If try_algo = 1 then \\code{W}, \\code{beta} and \\code{sigma2} are
      initialized by segmenting uniformly into \\code{K} contiguous segments
      the response Y. Otherwise, \\code{W}, \\code{beta} and \\code{sigma2} are
      initialized by segmenting randomly into \\code{K} segments the response Y."
      if (try_algo == 1) { # Uniform segmentation into K contiguous segments, and then a regression

        # Initialization of W
        W <<- zeros(q + 1, K - 1)

        zi <- round(fData$m / K) - 1

        beta <<- matrix(NA, p + 1, K)

        for (k in 1:K) {
          i <- (k - 1) * zi + 1
          j <- k * zi
          yij <- fData$Y[i:j]

          Phi_ij <- phi$XBeta[i:j, ]

          bk <-  solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% yij
          beta[, k] <<- bk

          if (variance_type == variance_types$homoskedastic) {
            sigma2 <<- matrix(1)
          }
          else{
            sigma2[k] <<- var(yij)
          }
        }
      } else{# Random segmentation into K contiguous segments, and then a regression

        # Initialization of W
        W <<- rand(q + 1, K - 1)

        Lmin <- round(fData$m / (K + 1)) # Minimum number of points in a segment
        tk_init <- zeros(K, 1)
        if (K == 1) {
          tk_init[2] = fData$m
        } else {
          K_1 <- K
          for (k in 2:K) {
            K_1 <- (K_1 - 1)
            temp <-
              (tk_init[k - 1] + Lmin):(fData$m - (K_1 * Lmin))

            ind <- sample(length(temp))

            tk_init[k] <- temp[ind[1]]
          }
          tk_init[K + 1] <- fData$m
        }

        beta <<- matrix(NA, p + 1, K)
        for (k in 1:K) {
          i <- tk_init[k] + 1
          j <- tk_init[k + 1]
          yij <- fData$Y[i:j]
          Phi_ij <- phi$XBeta[i:j, ]
          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% yij
          beta[, k] <<- bk

          if (variance_type == variance_types$homoskedastic) {
            sigma2 <<- var(fData$Y)
          }
          else{
            sigma2[k] <<- 1
          }
        }
      }
    },

    MStep = function(statRHLP, verbose_IRLS) {
      "Method used in the EM algorithm to learn the parameters of the RHLP model
      based on statistics provided by \\code{statRHLP}."
      # Maximization w.r.t betak and sigmak (the variances)
      if (variance_type == variance_types$homoskedastic) {
        s = 0
      }
      for (k in 1:K) {
        weights <- statRHLP$tau_ik[, k] # Post prob of each component k (dimension nx1)
        nk <- sum(weights) # Expected cardinal number of class k

        Xk <- phi$XBeta * (sqrt(weights) %*% ones(1, p + 1)) # [m*(p+1)]
        yk <- fData$Y * (sqrt(weights)) # Dimension: (nx1).*(nx1) = (nx1)

        M <- t(Xk) %*% Xk
        epps <- 1e-9
        M <- M + epps * diag(p + 1)

        beta[, k] <<- solve(M) %*% t(Xk) %*% yk # Maximization w.r.t betak
        z <- sqrt(weights) * (fData$Y - phi$XBeta %*% beta[, k])

        # Maximisation w.r.t sigmak (the variances)
        sk <- t(z) %*% z

        if (variance_type == variance_types$homoskedastic) {
          s <- s + sk
          sigma2 <<- s / fData$m
        } else {
          sigma2[k] <<- sk / nk
        }
      }

      # Maximization w.r.t W
      #  IRLS : Iteratively Reweighted Least Squares (for IRLS, see the IJCNN 2009 paper)
      res_irls <- IRLS(phi$Xw, statRHLP$tau_ik, ones(nrow(statRHLP$tau_ik), 1), W, verbose_IRLS)

      W <<- res_irls$W
      pi_ik <- res_irls$piik
      reg_irls <- res_irls$reg_irls
    }
  )
)
