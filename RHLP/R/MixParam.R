source("R/enums.R")
source("R/utils.R")

MixParam <- setRefClass(
  "MixParam",
  fields = list(
    Wk = "matrix",
    betak = "matrix",
    sigmak = "matrix"
  ),
  methods = list(
    initParam = function(mix, phi, mixOptions, try_algo = 1){
      m <- length(mix$y)

      if (try_algo==1){
        #Initialization of W
        Wk <<- zeros(mix$q+1,mix$K-1)

        # decoupage de l'echantillon (signal) en K segments
        zi <- round(m/mix$K)-1

        betak <<- matrix(NA, mix$p+1, mix$K)

        for (k in 1:mix$K){
          i <- (k-1)*zi+1
          j <- k*zi
          yij <- mix$y[i:j]
          yij <- matrix(t(yij), ncol = 1)
          Phi_ij <- phi$phiBeta[i:j,]

          bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%yij
          betak[,k] <<- bk

          if (mixOptions$variance_type == variance_types$homoskedastic){
            sigmak <<- matrix(1)
          }
          else{
            sigmak[k] <<- var(yij)
          }
        }
      }
      else{ # initialisation aléatoire
        #Initialization of W
        Wk <<- rand(mix$q+1,mix$K-1)

        Lmin <- round(m/(mix$K+1)) #nbr pts min dans un segments
        tk_init <- zeros(1,mix$K+1)
        K_1 <- mix$K
        for (k in 2:mix$K) {
          K_1 <- K_1-1;
          temp <- (tk_init[k-1] + Lmin) : (m - (K_1*Lmin))
          ind <- sample(length(temp));
          tk_init[k] <- temp[ind[1]]
        }
        tk_init[mix$K+1] <- m

        beta_k <- matrix(NA, mix$p+1, mix$K)

        for (k in 1:mix$K){
          i <- tk_init[k] + 1
          j <- tk_init[k+1]
          yij <- mix$y[i:j]
          yij <- matrix(t(yij), ncol = 1)
          Phi_ij <- phi$phiBeta[i:j,]

          bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%yij
          betak[,k] <<- bk

          if (mixOptions$variance_type == variance_types$homoskedastic){
            sigmak <<- var(mix$y)
          }
          else{
            sigmak[k] <<- 1
          }
        }
      }
    },

    MStep = function(mixModel, mixStats, phi, mixOptions){
      # M-Step
      # Maximization w.r.t betak and sigmak (the variances)
      if (mixOptions$variance_type == variance_types$homoskedastic) {s = 0}
      for (k in 1:mixModel$K){
        weights <- mixStats$tik[,k]; # post prob of each component k (dimension nx1)
        nk <- sum(weights); # expected cardinal numnber of class k

        Xk <- phi$phiBeta*(sqrt(weights)%*%ones(1,mixModel$p+1)); #[m*(p+1)]
        yk <- mixModel$y*(sqrt(weights));# dimension :(nx1).*(nx1) = (nx1)

        M <- t(Xk)%*%Xk
        epps <- 1e-9
        M <- M+epps*diag(mixModel$p+1);
        betak[,k] <<- solve(M)%*%t(Xk)%*%yk # Maximization w.r.t betak
        z <- sqrt(weights)*(mixModel$y-phi$phiBeta%*%betak[,k])
        # Maximisation w.r.t sigmak (the variances)
        priorsigma =  0; #1e-5;
        if (mixOptions$variance_type == variance_types$homoskedastic){
          sk <- t(z)%*%z
          s <- s+sk;
          sigmak <<- s/mixModel$m;
        }else{
          sigmak[k] <<- t(z)%*%z/nk  + priorsigma
        }
      }

      # Maximization w.r.t W
      # ----------------------------------%
      #  IRLS : Iteratively Reweighted Least Squares (for IRLS, see the IJCNN 2009 paper)
      res_irls <- IRLS_MixFRHLP(mixStats$tik, phi$phiW, Wk, verbose_IRLS = mixOptions$verbose_IRLS, piik_len = mixModel$m)

      Wk <<- res_irls[[1]]
      piik <- res_irls[[2]]
      reg_irls <- res_irls[[3]]
    }
  )
)

MixParam<-function(mixModel, options){
  #mixModel <- mixModel
  Wk <- matrix(0,mixModel$p+1, mixModel$K-1)
  betak <- matrix(NA, mixModel$p+1, mixModel$K)
  if (options$variance_type == variance_types$homoskedastic){
    # todo: verify the dimensions
    sigmak <- matrix(NA, mixModel$K)
  }
  else{
    sigmak <- matrix(NA, mixModel$K)
  }
  new("MixParam", Wk=Wk, betak=betak, sigmak=sigmak)#, mixModel = mixModel)
}
