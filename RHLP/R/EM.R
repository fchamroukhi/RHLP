source("R/IRLS.R")
source("R/StatRHLP.R")
source("R/model_logit.R")
source("R/utils.R")
source("R/regressionDesigner.R")

initParam <- function(modelRHLP, fData, phi, try_algo = 1){
  betak <- matrix(NA, modelRHLP$p+1, modelRHLP$K)
  if (modelRHLP$variance_type == variance_types$homoskedastic){
    sigmak <- matrix(1)
  }
  else{
    sigmak <- matrix(NA, K)
  }

  if (try_algo == 1){
    #Initialization of W
    Wk <- zeros(modelRHLP$q+1,modelRHLP$K-1)

    # decoupage de l'echantillon (signal) en K segments
    zi <- round(fData$m/modelRHLP$K) - 1

    for (k in 1:modelRHLP$K){
      i <- (k-1)*zi+1
      j <- k*zi
      yij <- fData$Y[i:j]

      Phi_ij <- phi$XBeta[i:j,]

      bk <-  solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%yij
      betak[, k] <- bk

      if (modelRHLP$variance_type == variance_types$homoskedastic){
        sigmak <- matrix(1)
      }
      else{
        sigmak[k] <- var(yij)
      }
    }
  }
  else{ # initialisation aléatoire
    #Initialization of W
    Wk <- rand(modelRHLP$q+1,modelRHLP$K-1)

    Lmin <- round(fData$m/(modelRHLP$K+1)) #nbr pts min dans un segments
    tk_init <- zeros(modelRHLP$K, 1)
    if (modelRHLP$K==1){
      tk_init[2] = m
    }
    else{
      K_1 <- modelRHLP$K
      for (k in 2:modelRHLP$K) {
        K_1 <- (K_1-1)
        temp <- (tk_init[k-1] + Lmin) : (fData$m - (K_1*Lmin))

        ind <- sample(length(temp));
        tk_init[k] <- temp[ind[1]]
      }
      tk_init[modelRHLP$K+1] <- fData$m
    }


    for (k in 1:modelRHLP$K){
      i <- tk_init[k] + 1
      j <- tk_init[k+1]
      yij <- fData$Y[i:j]
      #yij <- matrix(t(yij), ncol = 1)
      Phi_ij <- phi$XBeta[i:j,]
      bk <- solve(t(Phi_ij)%*%Phi_ij)%*%t(Phi_ij)%*%yij
      betak[,k] <- bk

      if (modelRHLP$variance_type == variance_types$homoskedastic){
        sigmak <- var(fData$Y)
      }
      else{
        sigmak[k] <- 1
      }
    }
  }

  return(list(W = Wk, beta = betak, sigma = sigmak))
}






EM <- function(K, p, q, variance_type, fData, n_tries, max_iter, threshold, verbose, verbose_IRLS){
  modelRHLP = ModelRHLP(K, p, q, variance_type)

  phi <- designmatrix(fData$X, modelRHLP$p, modelRHLP$q)

  top <- 0
  try_EM <- 0
  best_loglik <- -Inf
  cpu_time_all <- c()

  while(try_EM < n_tries){
    try_EM <- try_EM+1
    message("EM try nr ",try_EM)
    time <- Sys.time()

    #############################################
    # Initializations
    #############################################
    theta <- initParam(modelRHLP, fData, phi,  try_EM)
    modelRHLP$W <- theta$W
    modelRHLP$beta <- theta$beta
    modelRHLP$sigma <- theta$sigma

    iter <- 0
    converge <- FALSE
    prev_loglik <- -Inf

    stat <- StatRHLP(modelRHLP, fData$m)

    while(!converge && (iter <= max_iter)){
      #############################################
      # E Step
      #############################################
      stat$piik <- modele_logit(modelRHLP$W, phi$Xw)[[1]]
      for (k in (1:modelRHLP$K)){
        muk <- phi$XBeta%*%modelRHLP$beta[,k]
        if (modelRHLP$variance_type == variance_types$homoskedastic){
          sigmak <-  modelRHLP$sigma
        }
        else{
          sigmak <- modelRHLP$sigma[k]
        }
        z <- ((fData$Y-muk)^2)/sigmak
        stat$log_piik_fik[,k] <- log(stat$piik[,k]) - (0.5 * ones(fData$m,1) * (log(2*pi) + log(sigmak))) - (0.5*z)

      }

      stat$log_piik_fik <- pmax(stat$log_piik_fik, log(.Machine$double.xmin))
      piik_fik <- exp(stat$log_piik_fik)
      fxi <- rowSums(piik_fik)
      log_fxi <- log(fxi)
      stat$log_sum_piik_fik <- matrix(log(rowSums(piik_fik)))
      log_tik <- stat$log_piik_fik - stat$log_sum_piik_fik%*%ones(1,modelRHLP$K)
      stat$tik <- normalize(exp(log_tik),2)$M






      #############################################
      # M Step
      #############################################
      # Maximization w.r.t betak and sigmak (the variances)
      if (modelRHLP$variance_type == variance_types$homoskedastic) {s = 0}
      for (k in 1:modelRHLP$K){
        weights <- stat$tik[,k]; # post prob of each component k (dimension nx1)
        nk <- sum(weights); # expected cardinal numnber of class k

        Xk <- phi$XBeta*(sqrt(weights)%*%ones(1,modelRHLP$p+1)); #[m*(p+1)]
        yk <- fData$Y*(sqrt(weights));# dimension :(nx1).*(nx1) = (nx1)

        M <- t(Xk)%*%Xk
        epps <- 1e-9
        M <- M+epps*diag(modelRHLP$p+1);
        modelRHLP$beta[,k] <- solve(M)%*%t(Xk)%*%yk # Maximization w.r.t betak
        z <- sqrt(weights)*(fData$Y-phi$XBeta%*%modelRHLP$beta[,k])
        # Maximisation w.r.t sigmak (the variances)
        priorsigma =  0; #1e-5;
        if (modelRHLP$variance_type == variance_types$homoskedastic){
          sk <- t(z)%*%z
          s <- s+sk;
          modelRHLP$sigma <- s/modelRHLP$m;
        }else{
          modelRHLP$sigma[k] <- t(z)%*%z/nk  + priorsigma
        }
      }

      # Maximization w.r.t W
      # ----------------------------------%
      #  IRLS : Iteratively Reweighted Least Squares (for IRLS, see the IJCNN 2009 paper)
      res_irls <- IRLS(stat$tik, phi$Xw, modelRHLP$W, verbose_IRLS = verbose_IRLS, piik_len = mixModel$m)

      modelRHLP$W <- res_irls$W
      reg_irls <- res_irls$reg_irls



      stat$computeLikelihood(reg_irls)



      # FIN EM

      iter <- iter + 1
      if (verbose){
        message("EM     : Iteration : ",iter, "  log-likelihood : "  , stat$log_lik)
      }
      if (prev_loglik - stat$log_lik > 1e-5){
        message("!!!!! EM log-likelihood is decreasing from ", prev_loglik, "to ",  stat$log_lik)
        top <- top + 1
        if (top > 20) break
      }

      # TEST OF CONVERGENCE
      converge <- abs((stat$log_lik - prev_loglik)/prev_loglik) <= threshold
      if (is.na(converge)) {converge <- FALSE} # basicly for the first iteration when prev_loglik is Inf

      prev_loglik <- stat$log_lik
      stat$stored_loglik[iter] <- stat$log_lik
    }# FIN EM LOOP

    cpu_time_all[try_EM] <- Sys.time()-time

    # at this point we have computed param and stat that contains all the information

    if (stat$log_lik > best_loglik){
      statSolution <- stat$copy()
      paramSolution <- modelRHLP$copy()
      if (modelRHLP$K==1){
        statSolution$tik <- matrix(stat$tik, nrow = fData$m, ncol = 1)
        statSolution$piik <- matrix(stat$piik, nrow = fData$m, ncol = 1)
      }
      else{
        statSolution$tik <- stat$tik[1:fData$m,]
        statSolution$piik <- stat$piik[1:fData$m,]
      }

      best_loglik <- stat$log_lik
    }
    if (n_tries > 1){
      message("max value: ", stat$log_lik)
    }
  }

  # Computation of c_ig the hard partition of the curves and klas
  statSolution$MAP()

  if (n_tries > 1){
    message("max value: ", statSolution$log_lik)
  }


  # FINISH computation of statSolution
  statSolution$computeStats(modelRHLP, fData$m, phi, cpu_time_all)

  return(list(paramSolution=paramSolution, statSolution=statSolution))
}
