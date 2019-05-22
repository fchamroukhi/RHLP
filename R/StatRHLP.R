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
