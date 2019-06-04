#' A Reference Class which contains statistics of a RHLP model.
#'
#' StatRHLP contains all the parameters of a [RHLP][ParamRHLP] model.
#'
#' @field pi_ik Matrix of size \eqn{(m, K)} representing the probabilities
#' \eqn{P(zi = k; W) = P(z_{ik} = 1; W)}{P(zi = k; W) = P(z_ik = 1; W)} of the
#' latent variable \eqn{zi,\ i = 1,\dots,m}{zi, i = 1,\dots,m}.
#' @field z_ik Hard segmentation logical matrix of dimension \eqn{(m, K)}
#' obtained by the Maximum a posteriori (MAP) rule:
#' \eqn{z_{ik} = 1 \ \textrm{if} \ z_{ik} = \textrm{arg} \ \textrm{max}_{k} \
#' P(z_i = k | Y, W, \beta) = tau_ik;\ 0 \ \textrm{otherwise}}{z_ik = 1 if z_ik =
#' arg max_k P(z_i = k | Y, W, \beta) = tau_ik; 0 otherwise}, \eqn{k = 1,\dots,K}.
#' @field klas Column matrix of the labels issued from `z_ik`. Its elements are
#' \eqn{klas(i) = k}, \eqn{k = 1,\dots,K}.
#' @field Ex Column matrix of dimension \emph{m}. `Ex` is the curve expectation
#' : sum of the polynomial components \eqn{\beta_{k} \times X_{i}}{\betak x X_i}
#' weighted by the logistic probabilities `pi_ik`:
#' \eqn{Ey(i) = \sum_{k = 1}^{K} pi_ik \times \beta_{k} \times X_{i}}{Ey(i) =
#' \sum_{k=1}^K pi_ik x \betak x X_i}, \eqn{i = 1,\dots,m}.
#' @field loglik Numeric. Log-likelihood of the RHLP model.
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
#' criterion. The formula is \eqn{AIC = log\_lik - nu}{AIC = loglik - nu}.
#' @field cpu_time Numeric. Average executon time of a EM step.
#' @field log_piik_fik Matrix of size \eqn{(m, K)} giving the values of the
#' logarithm of the joint probability
#' \eqn{P(Y_{i}, \ zi = k)}{P(Yi, zi = k)}, \eqn{i = 1,\dots,m}.
#' @field log_sum_piik_fik Column matrix of size \emph{m} giving the values of
#' \eqn{\sum_{k = 1}^{K} \textrm{log} P(Y_{i}, \ zi = k)}{\sum_{k = 1}^{K} log
#' P(Yi, zi = k)}, \eqn{i = 1,\dots,m}.
#' @field tau_ik Matrix of size \eqn{(m, K)} giving the posterior probability that
#' \eqn{Y_{i}}{Yi} originates from the \eqn{k}-th regression model
#' \eqn{P(zi = k | Y, W, \beta)}.
#' @field polynomials Matrix of size \eqn{(m, K)} giving the values of
#' \eqn{\beta_{k} \times X_{i}}{\betak x X_i}, \eqn{i = 1,\dots,m}.
#' @seealso [ParamRHLP], [FData]
#' @export
StatRHLP <- setRefClass(
  "StatRHLP",
  fields = list(
    pi_ik = "matrix",
    z_ik = "matrix",
    klas = "matrix",
    Ex = "matrix",
    loglik = "numeric",
    com_loglik = "numeric",
    stored_loglik = "list",
    stored_com_loglik = "list",
    BIC = "numeric",
    ICL = "numeric",
    AIC = "numeric",
    cpu_time = "numeric",
    log_piik_fik = "matrix",
    log_sum_piik_fik = "matrix",
    tau_ik = "matrix",
    polynomials = "matrix"
  ),
  methods = list(

    initialize = function(paramRHLP = ParamRHLP(fData = FData(numeric(1), matrix(1)), K = 1, p = 2, q = 1, variance_type = 1)) {

      pi_ik <<- matrix(NA, paramRHLP$fData$m, paramRHLP$K)
      z_ik <<- matrix(NA, paramRHLP$fData$m, paramRHLP$K)
      klas <<- matrix(NA, paramRHLP$fData$m, 1)
      Ex <<- matrix(NA, paramRHLP$fData$m, 1)
      loglik <<- -Inf
      com_loglik <<- -Inf
      stored_loglik <<- list()
      stored_com_loglik <<- list()
      BIC <<- -Inf
      ICL <<- -Inf
      AIC <<- -Inf
      cpu_time <<- Inf
      log_piik_fik <<- matrix(0, paramRHLP$fData$m, paramRHLP$K)
      log_sum_piik_fik <<- matrix(NA, paramRHLP$fData$m, 1)
      tau_ik <<- matrix(0, paramRHLP$fData$m, paramRHLP$K)
      polynomials <<- matrix(NA, paramRHLP$fData$m, paramRHLP$K)

    },

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
      N <- nrow(pi_ik)
      K <- ncol(pi_ik)
      ikmax <- max.col(pi_ik)
      ikmax <- matrix(ikmax, ncol = 1)
      z_ik <<- ikmax %*% ones(1, K) == ones(N, 1) %*% (1:K)
      klas <<- ones(N, 1)
      for (k in 1:K) {
        klas[z_ik[, k] == 1] <<- k
      }
    },

    computeLikelihood = function(reg_irls) {
      loglik <<- sum(log_sum_piik_fik) + reg_irls

    },

    computeStats = function(paramRHLP, cpu_time_all) {

      polynomials <<- paramRHLP$phi$XBeta %*% paramRHLP$beta
      Ex <<- matrix(rowSums(pi_ik * polynomials))

      cpu_time <<- mean(cpu_time_all)

      BIC <<- loglik - (paramRHLP$nu * log(paramRHLP$fData$m) / 2)
      AIC <<- loglik - paramRHLP$nu

      zik_log_alphag_fg_xij <- (z_ik) * (log_piik_fik)
      com_loglik <<- sum(rowSums(zik_log_alphag_fg_xij))
      ICL <<- com_loglik - paramRHLP$nu * log(paramRHLP$fData$m) / 2

    },

    EStep = function(paramRHLP) {
      "Method used in the EM algorithm to update statistics based on parameters
      provided by \\code{paramRHLP} (prior and posterior probabilities)."
      pi_ik <<- multinomialLogit(paramRHLP$W, paramRHLP$phi$Xw, ones(paramRHLP$fData$m, paramRHLP$K), ones(paramRHLP$fData$m, 1))$piik

      for (k in (1:paramRHLP$K)) {
        muk <- paramRHLP$phi$XBeta %*% paramRHLP$beta[, k]

        if (paramRHLP$variance_type == variance_types$homoskedastic) {
          sigmak <-  paramRHLP$sigma[1]
        } else {
          sigmak <- paramRHLP$sigma[k]
        }
        z <- ((paramRHLP$fData$Y - muk) ^ 2) / sigmak
        log_piik_fik[, k] <<-
          log(pi_ik[, k]) - (0.5 * ones(paramRHLP$fData$m, 1) * (log(2 * pi) + log(sigmak))) - (0.5 * z)
      }

      log_piik_fik <<- pmax(log_piik_fik, log(.Machine$double.xmin))
      piik_fik <- exp(log_piik_fik)
      fxi <- rowSums(piik_fik)
      log_fxi <- log(fxi)
      log_sum_piik_fik <<- matrix(log(rowSums(piik_fik)))
      log_tik <- log_piik_fik - log_sum_piik_fik %*% ones(1, paramRHLP$K)
      tau_ik <<- normalize(exp(log_tik), 2)$M
    }
  )
)
