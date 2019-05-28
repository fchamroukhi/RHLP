#' A Reference Class which contains statistics of a RHLP model.
#'
#' StatRHLP contains all the parameters of a [RHLP][ModelRHLP] model.
#'
#' @usage NULL
#' @field piik Matrix of size \eqn{(m, K)} representing the probabilities
#' \eqn{P(zi = k; W) = P(z_{ik} = 1; W)}{P(zi = k; W) = P(z_ik = 1; W)} of the
#' latent variable \eqn{zi,\ i = 1,\dots,m}{zi, i = 1,\dots,m}.
#' @field z_ik Hard segmentation logical matrix of dimension \eqn{(m, K)}
#' obtained by the Maximum a posteriori (MAP) rule:
#' \eqn{z_{ik} = 1 \ \textrm{if} \ z_{ik} = \textrm{arg} \ \textrm{max}_{k} \
#' P(z_i = k | Y, W, \beta) = tik;\ 0 \ \textrm{otherwise}}{z_ik = 1 if z_ik =
#' arg max_k P(z_i = k | Y, W, \beta) = tik; 0 otherwise}, \eqn{k = 1,\dots,K}.
#' @field klas Column matrix of the labels issued from `z_ik`. Its elements are
#' \eqn{klas(i) = k}, \eqn{k = 1,\dots,K}.
#' @field Ex Column matrix of dimension \emph{m}. `Ex` is the curve expectation
#' : sum of the polynomial components \eqn{\beta_{k} \times X_{i}}{\betak x X_i}
#' weighted by the logistic probabilities `piik`:
#' \eqn{Ey(i) = \sum_{k = 1}^{K} piik \times \beta_{k} \times X_{i}}{Ey(i) =
#' \sum_{k=1}^K piik x \betak x X_i}, \eqn{i = 1,\dots,m}.
#' @field log_lik Numeric. Log-likelihood of the RHLP model.
#' @field com_loglik Numeric. Complete log-likelihood of the RHLP model.
#' @field stored_loglik List. Stored values of the log-likelihood at each EM
#' iteration.
#' @field BIC Numeric. Value of the BIC (Bayesian Information Criterion)
#' criterion. The formula is \eqn{BIC = log\_lik - nu \times
#' \textrm{log}(m) / 2}{BIC = loglik - nu x log(m) / 2} with \emph{nu} the
#' degree of freedom of the RHLP model.
#' @field ICL Numeric. Value of the ICL (Integrated Completed Likelihood)
#' criterion. The formula is \eqn{ICL = com\_loglik - nu \times
#' \textrm{log}(m) / 2}{ICL = com_loglik - nu x log(m) / 2} with \emph{nu} the
#' degree of freedom of the RHLP model.
#' @field AIC Numeric. Value of the AIC (Akaike Information Criterion)
#' criterion. The formula is \eqn{AIC = log\_lik - nu}{AIC = log_lik - nu}.
#' @field cpu_time Numeric. Average executon time of a EM step.
#' @field log_piik_fik Matrix of size \eqn{(m, K)} giving the values of the
#' logarithm of the joint probability
#' \eqn{P(Y_{i}, \ zi = k)}{P(Yi, zi = k)}, \eqn{i = 1,\dots,m}.
#' @field log_sum_piik_fik Column matrix of size \emph{m} giving the values of
#' \eqn{\sum_{k = 1}^{K} \textrm{log} P(Y_{i}, \ zi = k)}{\sum_{k = 1}^{K} log
#' P(Yi, zi = k)}, \eqn{i = 1,\dots,m}.
#' @field tik Matrix of size \eqn{(m, K)} giving the posterior probability that
#' \eqn{Y_{i}}{Yi} originates from the \eqn{k}-th regression model
#' \eqn{P(zi = k | Y, W, \beta)}.
#' @field polynomials Matrix of size \eqn{(m, K)} giving the values of
#' \eqn{\beta_{k} \times X_{i}}{\betak x X_i}, \eqn{i = 1,\dots,m}.
#' @field weighted_polynomials Matrix of size \emph{(m, K)} giving the values
#' of \eqn{piik \times \beta_{k} \times X_{i}}{piik x \betak x X_i},
#' \eqn{i = 1,\dots,m}.
#' @seealso [FData]
#' @export
StatRHLP <- setRefClass(
  "StatRHLP",
  fields = list(
    piik = "matrix",
    z_ik = "matrix",
    klas = "matrix",
    Ex = "matrix",
    log_lik = "numeric",
    com_loglik = "numeric",
    stored_loglik = "list",
    stored_com_loglik = "list",
    BIC = "numeric",
    ICL = "numeric",
    AIC = "numeric",
    cpu_time = "numeric",
    log_piik_fik = "matrix",
    log_sum_piik_fik = "matrix",
    tik = "matrix",
    polynomials = "matrix",
    weighted_polynomials = "matrix"
  ),
  methods = list(
    MAP = function() {
      ###########################################################################
      # function [klas, Z] = MAP(PostProbs)
      #
      # calculate a partition by applying the Maximum A Posteriori Bayes
      # allocation rule
      #
      #
      # Inputs :
      #   PostProbs, a matrix of dimensions [n x K] of the posterior
      #  probabilities of a given sample of n observations arizing from K groups
      #
      # Outputs:
      #   klas: a vector of n class labels (z_1, ...z_n) where z_i =k \in {1,...K}
      #       klas(i) = arg   max (PostProbs(i,k)) , for all i=1,...,n
      #                     1<=k<=K
      #               = arg   max  p(zi=k|xi;theta)
      #                     1<=k<=K
      #               = arg   max  p(zi=k;theta)p(xi|zi=k;theta)/sum{l=1}^{K}p(zi=l;theta) p(xi|zi=l;theta)
      #                     1<=k<=K
      #
      #
      #       Z : Hard partition data matrix [nxK] with binary elements Zik such
      #       that z_ik =1 iff z_i = k
      #
      N <- nrow(piik)
      K <- ncol(piik)
      ikmax <- max.col(piik)
      ikmax <- matrix(ikmax, ncol = 1)
      z_ik <<- ikmax %*% ones(1, K) == ones(N, 1) %*% (1:K)
      klas <<- ones(N, 1)
      for (k in 1:K) {
        klas[z_ik[, k] == 1] <<- k
      }
    },

    computeLikelihood = function(reg_irls) {
      log_lik <<- sum(log_sum_piik_fik) + reg_irls

    },

    computeStats = function(modelRHLP, paramRHLP, phi, cpu_time_all) {
      polynomials <<- phi$XBeta %*% paramRHLP$beta
      weighted_polynomials <<- piik * polynomials
      Ex <<- matrix(rowSums(weighted_polynomials))

      cpu_time <<- mean(cpu_time_all)

      BIC <<- log_lik - (modelRHLP$nu * log(modelRHLP$m) / 2)
      AIC <<- log_lik - modelRHLP$nu

      zik_log_alphag_fg_xij <- (z_ik) * (log_piik_fik)
      com_loglik <<- sum(rowSums(zik_log_alphag_fg_xij))
      ICL <<- com_loglik - modelRHLP$nu * log(modelRHLP$m) / 2

    },

    EStep = function(modelRHLP, paramRHLP, phi) {

      piik <<- multinomialLogit(paramRHLP$W, phi$Xw, ones(modelRHLP$m, modelRHLP$K), ones(modelRHLP$m, 1))$piik

      for (k in (1:K)) {
        muk <- phi$XBeta %*% paramRHLP$beta[, k]

        if (modelRHLP$variance_type == variance_types$homoskedastic) {
          sigmak <-  paramRHLP$sigma[1]
        } else {
          sigmak <- paramRHLP$sigma[k]
        }
        z <- ((modelRHLP$Y - muk) ^ 2) / sigmak
        log_piik_fik[, k] <<-
          log(piik[, k]) - (0.5 * ones(modelRHLP$m, 1) * (log(2 * pi) + log(sigmak))) - (0.5 * z)
      }

      log_piik_fik <<- pmax(log_piik_fik, log(.Machine$double.xmin))
      piik_fik <- exp(log_piik_fik)
      fxi <- rowSums(piik_fik)
      log_fxi <- log(fxi)
      log_sum_piik_fik <<- matrix(log(rowSums(piik_fik)))
      log_tik <- log_piik_fik - log_sum_piik_fik %*% ones(1, modelRHLP$K)
      tik <<- normalize(exp(log_tik), 2)$M
    }
  )
)

StatRHLP <- function(modelRHLP) {
  piik <- matrix(NA, modelRHLP$m, modelRHLP$K)
  z_ik <- matrix(NA, modelRHLP$m, modelRHLP$K)
  klas <- matrix(NA, modelRHLP$m, 1)
  Ex <- matrix(NA, modelRHLP$m, 1)
  log_lik <- -Inf
  com_loglik <- -Inf
  stored_loglik <- list()
  stored_com_loglik <- list()
  BIC <- -Inf
  ICL <- -Inf
  AIC <- -Inf
  cpu_time <- Inf
  log_piik_fik <- matrix(0, modelRHLP$m, modelRHLP$K)
  log_sum_piik_fik <- matrix(NA, modelRHLP$m, 1)
  tik <- matrix(0, modelRHLP$m, modelRHLP$K)
  polynomials <- matrix(NA, modelRHLP$m, modelRHLP$K)
  weighted_polynomials <- matrix(NA, modelRHLP$m, modelRHLP$K)

  new(
    "StatRHLP",
    piik = piik,
    z_ik = z_ik,
    klas = klas,
    Ex = Ex,
    log_lik = log_lik,
    com_loglik = com_loglik,
    stored_loglik = stored_loglik,
    stored_com_loglik = stored_com_loglik,
    BIC = BIC,
    ICL = ICL,
    AIC = AIC,
    cpu_time = cpu_time,
    log_piik_fik = log_piik_fik,
    log_sum_piik_fik = log_sum_piik_fik,
    tik = tik,
    polynomials = polynomials,
    weighted_polynomials = weighted_polynomials
  )
}
