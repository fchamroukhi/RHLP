MixStats <- setRefClass(
  "MixStats",
  fields = list(
    h_ig = "matrix", # post probabilities            piik in the matlab code of RHLP
    z_ik = "matrix", # c_ig in MIX_RHLP
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
    log_sum_piik_fik = "matrix",
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
      z_ik <<- ikmax%*%ones(1,K) == ones(N,1)%*%(1:K) # partition_MAP
      klas <<- ones(N,1)
      for (k in 1:K){
        klas[z_ik[,k]==1] <<- k
      }
    },
    #######
    # compute loglikelihood
    #######
    computeLikelihood = function(reg_irls){
      log_lik <<- sum(log_sum_piik_fik) + reg_irls;
    },
    #######
    #
    #######
    #######
    # compute the final solution stats
    #######
    computeStats = function(mixModel, mixParam, phi, cpu_time_all){
      polynomials <<- phi$phiBeta %*% mixParam$betak
      weighted_polynomials <<- h_ig * polynomials
      Ex <<- matrix(rowSums(weighted_polynomials))

      cpu_time <<- mean(cpu_time_all)
      Psi <- c(as.vector(mixParam$Wk), as.vector(mixParam$betak), as.vector(mixParam$sigmak))
      BIC <<- log_lik - (nu*log(mixModel$m)/2)
      AIC <<- log_lik - nu


      zik_log_alphag_fg_xij <- (z_ik)*(log_piik_fik);
      com_loglik <<- sum(rowSums(zik_log_alphag_fg_xij));

      ICL <<- com_loglik - nu*log(mixModel$m)/2;

    },
    #######
    # EStep
    #######
    EStep = function(mixModel, mixParam, phi, variance_type){
      h_ig <<- modele_logit(mixParam$Wk, phi$phiW)[[1]]
      for (k in (1:K)){
        muk <- phi$phiBeta%*%mixParam$betak[,k]
        if (variance_type == variance_types$homoskedastic){
          sigmak <-  mixParam$sigmak[1]
        }
        else{
          sigmak <- mixParam$sigmak[k]
        }
        z <- ((mixModel$y-muk)^2)/sigmak
        log_piik_fik[,k] <<- log(h_ig[,k]) - (0.5 * ones(mixModel$m,1) * (log(2*pi) + log(sigmak))) - (0.5*z)

      }

      log_piik_fik <<- pmax(log_piik_fik, log(.Machine$double.xmin))
      piik_fik <- exp(log_piik_fik)
      fxi <- rowSums(piik_fik)
      log_fxi <- log(fxi)
      log_sum_piik_fik <<- matrix(log(rowSums(piik_fik)))
      log_tik <- log_piik_fik - log_sum_piik_fik%*%ones(1,mixModel$K)
      tik <<- normalize(exp(log_tik),2)$M
    }
  )
)


MixStats<-function(mixModel, options){
  h_ig <- matrix(NA,mixModel$m, mixModel$K)
  z_ik <- matrix(NA,mixModel$m, mixModel$K)
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
  log_sum_piik_fik <- matrix(NA,mixModel$m, 1)
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
  new("MixStats", h_ig=h_ig, z_ik=z_ik, klas=klas, nu=nu, Ex=Ex, log_lik=log_lik, com_loglik=com_loglik, stored_loglik=stored_loglik, stored_com_loglik=stored_com_loglik, BIC=BIC, ICL=ICL, AIC=AIC, cpu_time=cpu_time,
      log_piik_fik=log_piik_fik, log_sum_piik_fik = log_sum_piik_fik, tik=tik, polynomials=polynomials, weighted_polynomials=weighted_polynomials)
}
