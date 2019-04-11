source("R/enums.R")
source("R/utils.R")
source("R/IRLS.R")

ParamRHLP <- setRefClass(
  "ParamRHLP",
  fields = list(W = "matrix",
                beta = "matrix",
                sigma = "matrix"),
  methods = list(
    initParam = function(modelRHLP, phi, try_algo = 1) {
      if (try_algo == 1) {
        #Initialization of W
        W <<- zeros(modelRHLP$q + 1, modelRHLP$K - 1)

        # decoupage de l'echantillon (signal) en K segments
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
      }
      else{
        # initialisation aléatoire
        #Initialization of W
        W <<- rand(modelRHLP$q + 1, modelRHLP$K - 1)

        Lmin <-
          round(modelRHLP$m / (modelRHLP$K + 1)) #nbr pts min dans un segments
        tk_init <- zeros(modelRHLP$K, 1)
        if (modelRHLP$K == 1) {
          tk_init[2] = modelRHLP$m
        }
        else{
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
          #yij <- matrix(t(yij), ncol = 1)
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
      # M-Step
      # Maximization w.r.t betak and sigmak (the variances)
      if (modelRHLP$variance_type == variance_types$homoskedastic) {
        s = 0
      }
      for (k in 1:modelRHLP$K) {
        weights <-
          statRHLP$tik[, k]
        # post prob of each component k (dimension nx1)
        nk <- sum(weights)
        # expected cardinal numnber of class k

        Xk <-
          phi$XBeta * (sqrt(weights) %*% ones(1, modelRHLP$p + 1))
        #[m*(p+1)]
        yk <-
          modelRHLP$Y * (sqrt(weights))
        # dimension :(nx1).*(nx1) = (nx1)

        M <- t(Xk) %*% Xk
        epps <- 1e-9
        M <- M + epps * diag(modelRHLP$p + 1)

        beta[, k] <<-
          solve(M) %*% t(Xk) %*% yk # Maximization w.r.t betak
        z <- sqrt(weights) * (modelRHLP$Y - phi$XBeta %*% beta[, k])
        # Maximisation w.r.t sigmak (the variances)
        priorsigma =  0
        #1e-5;
        if (modelRHLP$variance_type == variance_types$homoskedastic) {
          sk <- t(z) %*% z
          s <- s + sk

          sigma <<- s / modelRHLP$m

        } else{
          sigma[k] <<- t(z) %*% z / nk  + priorsigma
        }
      }

      # Maximization w.r.t W
      # ----------------------------------%
      #  IRLS : Iteratively Reweighted Least Squares (for IRLS, see the IJCNN 2009 paper)
      res_irls <-
        IRLS(statRHLP$tik,
             phi$Xw,
             W,
             verbose_IRLS = verbose_IRLS,
             piik_len = modelRHLP$m)

      W <<- res_irls[[1]]
      piik <- res_irls[[2]]
      reg_irls <- res_irls[[3]]
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
  new("ParamRHLP", W = W, beta = beta, sigma = sigma)#, modelRHLP = modelRHLP)
}
