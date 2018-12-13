MixStats <- setRefClass(
  "MixStats",
  fields = list(
    h_ig = "matrix", # post probabilities            piik in the matlab code of RHLP
    klas = "matrix",
    # pi_jgk = "matrix",
    nu="numeric",
    Ex = "matrix",
    log_lik="numeric",
    com_loglik="numeric",
    stored_loglik = "list",
    stored_com_loglik = "list",
    BIC="numeric",
    ICL="numeric",
    AIC="numeric",
    cpu_time = "numeric",
    log_piik_fik="matrix",
    tik="matrix",
    polynomials="matrix",
    weighted_polynomials="matrix"
  ),
  methods=list(
    MAP = function(){
      "
      calcule une partition d'un echantillon par la regle du Maximum A Posteriori à partir des probabilites a posteriori
      Entrees : post_probas , Matrice de dimensions [n x K] des probabibiltes a posteriori (matrice de la partition floue)
           n : taille de l'echantillon
           K : nombres de classes
           klas(i) = arg   max (post_probas(i,k)) , for all i=1,...,n
                         1<=k<=K
                   = arg   max  p(zi=k|xi;theta)
                         1<=k<=K
                   = arg   max  p(zi=k;theta)p(xi|zi=k;theta)/sum{l=1}^{K}p(zi=l;theta) p(xi|zi=l;theta)
                         1<=k<=K
      Sorties : classes : vecteur collones contenant les classe (1:K)
           Z : Matrice de dimension [nxK] de la partition dure : ses elements sont zik, avec zik=1 si xi
           appartient à la classe k (au sens du MAP) et zero sinon.
      "
      N <- nrow(h_ig)
      K <- ncol(h_ig)
      ikmax <- max.col(h_ig)
      ikmax <- matrix(ikmax, ncol = 1)
      c_ig <- ikmax%*%ones(1,K) == ones(N,1)%*%(1:K) # partition_MAP
      klas <<- ones(N,1)
      for (k in 1:K){
        klas[c_ig[,k]==1] <<- k
      }
    },

    #######
    # EStep
    #######
    EStep = function(mixModel, mixParam, phi, variance_type){
      piik <- modele_logit(Winit,phi$phiW);
    }
  )
)


MixStats<-function(mixModel, options){
  h_ig <- matrix(NA,mixModel$m, mixModel$K)
  klas <- matrix(NA, mixModel$m, 1)
  Ex <- matrix(NA,mixModel$m, 1)
  log_lik <- -Inf
  com_loglik <- -Inf
  stored_loglik <- list()
  stored_com_loglik <- list()
  BIC <- -Inf
  ICL <- -Inf
  AIC <- -Inf
  cpu_time <- Inf
  log_piik_fik <- matrix(0, mixModel$m, mixModel$K)
  tik <- matrix(0, mixModel$m, mixModel$K)
  polynomials <- matrix(NA, mixModel$m, mixModel$K)
  weighted_polynomials <- matrix(NA, mixModel$m, mixModel$K)
  #number of free model parameters
  if (options$variance_type == variance_types$homoskedastic){
    nu <<- (mixModel$p+mixModel$q+3)*mixModel$K-(mixModel$q+1) - (mixModel$K-1)
  }
  else{
    nu <<- (mixModel$p+mixModel$q+3)*mixModel$K-(mixModel$q+1)
  }
  new("MixStats", h_ig=h_ig, klas=klas, nu=nu, Ex=Ex, log_lik=log_lik, com_loglik=com_loglik, stored_loglik=stored_loglik, stored_com_loglik=stored_com_loglik, BIC=BIC, ICL=ICL, AIC=AIC, cpu_time=cpu_time,
      log_piik_fik=log_piik_fik, tik=tik, polynomials=polynomials, weighted_polynomials=weighted_polynomials)
}
