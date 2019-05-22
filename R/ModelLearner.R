EM <- function(modelRHLP, n_tries = 1, max_iter = 1500, threshold = 1e-6, verbose = FALSE, verbose_IRLS = FALSE) {

  phi <- designmatrix(x = modelRHLP$X, p = modelRHLP$p, q = modelRHLP$q)

  top <- 0
  try_EM <- 0
  best_loglik <- -Inf
  cpu_time_all <- c()

  while (try_EM < n_tries) {
    try_EM <- try_EM + 1
    message("EM try nr ", try_EM)
    time <- Sys.time()

    # Initialization
    param <- ParamRHLP(modelRHLP)
    param$initParam(modelRHLP, phi, try_EM)
    iter <- 0
    converge <- FALSE
    prev_loglik <- -Inf

    stat <- StatRHLP(modelRHLP)

    while (!converge && (iter <= max_iter)) {
      stat$EStep(modelRHLP, param, phi)

      reg_irls <- param$MStep(modelRHLP, stat, phi, verbose_IRLS)
      stat$computeLikelihood(reg_irls)

      iter <- iter + 1
      if (verbose) {
        message("EM     : Iteration : ", iter, "  log-likelihood : "  , stat$log_lik)
      }
      if (prev_loglik - stat$log_lik > 1e-5) {
        message("!!!!! EM log-likelihood is decreasing from ", prev_loglik, "to ", stat$log_lik)
        top <- top + 1
        if (top > 20)
          break
      }

      # Test of convergence
      converge <- abs((stat$log_lik - prev_loglik) / prev_loglik) <= threshold
      if (is.na(converge)) {
        converge <- FALSE
      } # Basically for the first iteration when prev_loglik is Inf

      prev_loglik <- stat$log_lik
      stat$stored_loglik[iter] <- stat$log_lik
    } # End of the EM loop

    cpu_time_all[try_EM] <- Sys.time() - time

    if (stat$log_lik > best_loglik) {
      statSolution <- stat$copy()
      paramSolution <- param$copy()
      if (modelRHLP$K == 1) {
        statSolution$tik <- matrix(stat$tik, nrow = modelRHLP$m, ncol = 1)
        statSolution$piik <- matrix(stat$piik, nrow = modelRHLP$m, ncol = 1)
      }
      else{
        statSolution$tik <- stat$tik[1:modelRHLP$m, ]
        statSolution$piik <- stat$piik[1:modelRHLP$m, ]
      }

      best_loglik <- stat$log_lik
    }
    if (n_tries > 1) {
      message("max value: ", stat$log_lik)
    }
  }

  # Computation of c_ig the hard partition of the curves and klas
  statSolution$MAP()

  if (n_tries > 1) {
    message("max value: ", statSolution$log_lik)
  }

  # Finish the computation of statistics
  statSolution$computeStats(modelRHLP, paramSolution, phi, cpu_time_all)

  return(FittedRHLP(modelRHLP, paramSolution, statSolution))
}
