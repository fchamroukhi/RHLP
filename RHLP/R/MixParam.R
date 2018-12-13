source("R/enums.R")
source("R/utils.R")

MixParam <- setRefClass(
  "MixParam",
  fields = list(
    Wk = "array",
    betak = "array",
    sigmak = "matrix",
    piik = "array"
  ),
  methods = list(

    # y,K,XBeta,type_variance, try_EM
    initRegressionParam = function(Xg, g, K, p, phiBeta, variance_type, try_algo){
       n <- nrow(Xg)
       m <- ncol(Xg)
       if (try_algo==1){
          # decoupage de l'echantillon (signal) en K segments
          zi <- round(m/K)-1

          beta_k <- matrix(NA, p+1, K)
          sigma <- c()

          for (k in 1:K){
            i <- (k-1)*zi+1
            j <- k*zi
            Xij <- Xg[,i:j]
            Xij <- matrix(t(Xij), ncol = 1)
            phi_ij <- phiBeta[i:j,]
            Phi_ij <- repmat(phi_ij,n,1)

            bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%Xij
            beta_k[,k] <- bk

            if (variance_type == variance_types$common){
              sigma <- var(Xij)
            }
            else{
              mk <- j-i+1 #length(Xij);
              z <- Xij - Phi_ij %*% bk;
              sk <- t(z) %*% z/(n*mk);
              sigma[k] <- sk;
            }
          }
       }
       else{ # initialisation aléatoire
         Lmin <- round(m/(K+1)) #nbr pts min dans un segments
         tk_init <- zeros(1,K+1)
         K_1 <- K
         for (k in 2:K) {
           K_1 <- K_1-1;
           temp <- (tk_init[k-1] + Lmin) : (m - (K_1*Lmin))
           ind <- sample(length(temp));
           tk_init[k] <- temp[ind[1]]
         }
         tk_init[K+1] <- m

         beta_k <- matrix(NA, p+1, K)
         sigma <- c()
         for (k in 1:K){
           i <- tk_init[k] + 1
           j <- tk_init[k+1]
           Xij <- Xg[,i:j]
           Xij <- matrix(t(Xij), ncol = 1)
           phi_ij <- phiBeta[i:j,]
           Phi_ij <- repmat(phi_ij,n,1)


           bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%Xij
           beta_k[,k] <- bk

           if (variance_type == variance_types$common){
             sigma <- var(Xij)
           }
           else{
             mk <- j-i+1 #length(Xij);
             z <- Xij - Phi_ij %*% bk;
             sk <- t(z) %*% z/(n*mk);
             sigma[k] <- sk;
           }
         }
       }

       betag[,,g] <<- beta_k
       if (variance_type == variance_types$common){
         sigmag[g] <<- sigma
       }
       else{
         sigmag[,g] <<- sigma
       }
    }



  )
)

MixParam<-function(mixModel, options){
  #mixModel <- mixModel
  Wk <- array(0,dim=c(mixModel$q+1, mixModel$K-1, mixModel$G))
  betak <- array(NA, dim=c(mixModel$p+1, mixModel$K, mixModel$G))
  if (options$variance_type == variance_types$homoskedastic){
    sigmak <- matrix(NA, mixModel$G)
  }
  else{
    sigmak <- matrix(NA, mixModel$K, mixModel$G)
  }
  piik <- array(0, dim=c(mixModel$m*mixModel$n, mixModel$K, mixModel$G))
  new("MixParam", Wk=Wk, betak=betak, sigmak=sigmak, piik=piik)#, mixModel = mixModel)
}
