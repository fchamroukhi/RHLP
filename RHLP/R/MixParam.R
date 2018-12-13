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
        sigmak <<- c()
        for (k in 1:mix$K){
          i <- tk_init[k] + 1
          j <- tk_init[k+1]
          yij <- mix$y[,i:j]
          yij <- matrix(t(yij), ncol = 1)
          Phi_ij <- phi$phiBeta[i:j,]

          bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%yij
          betak[,k] <<- bk

          if (mixOptions$variance_type == variance_types$homoskedastic){
            sigmak <<- var(y)
          }
          else{
            sigmak[k] <<- matrix(1);
          }
        }
      }
    },

    MStep = function(mixModel, mixStats, phi, mixOptions){

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
