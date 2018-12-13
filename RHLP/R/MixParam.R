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

    #
    initParam = function(mix, phi, mixOptions, try_algo = 1){
      y<-mix$y
      K<-mix$K
      variance_type <-mixOptions$variance_type
      p<- mix$p
       m <- length(y)
       if (try_algo==1){
          # decoupage de l'echantillon (signal) en K segments
          zi <- round(m/K)-1

          betak <<- matrix(NA, p+1, K)


          for (k in 1:K){
            i <- (k-1)*zi+1
            j <- k*zi
            yij <- y[i:j]
            yij <- matrix(t(yij), ncol = 1)
            Phi_ij <- phi$phiBeta[i:j,]

            bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%yij
            betak[,k] <<- bk

            if (variance_type == variance_types$homoskedastic){
              sigmak <<- matrix(1)
            }
            else{
              sigmak[k] <<- var(yij)
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
         sigmak <<- c()
         for (k in 1:K){
           i <- tk_init[k] + 1
           j <- tk_init[k+1]
           yij <- y[,i:j]
           yij <- matrix(t(yij), ncol = 1)
           Phi_ij <- phi$phiBeta[i:j,]

           bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%yij
           betak[,k] <<- bk

           if (variance_type == variance_types$homoskedastic){
             sigmak <<- var(y)
           }
           else{
             sigmak[k] <<- matrix(1);
           }
         }
       }
    }



  )
)

MixParam<-function(mixModel, options){
  #mixModel <- mixModel
  Wk <- matrix(0,mixModel$p+1, mixModel$K)
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
