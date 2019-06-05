#' emRHLP is used to fit a RHLP model.
#'
#' emRHLP is used to fit a RHLP model. The estimation method is performed by
#' the Expectation-Maximization algorithm.
#'
#' @details emRHLP function is based on the EM algorithm. This functions starts
#' with an initialization of the parameters done by the method `initParam` of
#' the class [ParamRHLP][ParamRHLP], then it alternates between a E-Step
#' (method of the class [StatRHLP][StatRHLP]) and a M-Step (method of the class
#' [ParamRHLP][ParamRHLP]) until convergence (until the absolute difference of
#' log-likelihood between two steps of the EM algorithm is less than the
#' `threshold` parameter).
#'
#' @param X Numeric vector of length \emph{m} representing the covariates.
#' @param Y Matrix of size \eqn{(n, m)} representing \emph{n} functions of `X`
#' observed at points \eqn{1,\dots,m}.
#' @param K The number of regimes (mixture components).
#' @param p The order of the polynomial regression.
#' @param q The dimension of the logistic regression. For the purpose of
#' segmentation, it must be set to 1.
#' @param variance_type Numeric indicating if the model is homoskedastic
#' (`variance_type` = 1) or heteroskedastic (`variance_type` = 2).
#' @param n_tries Number of times EM algorithm will be launched.
#' The solution providing the highest log-likelihood will be returned.
#'
#' If `n_tries` > 1, then for the first pass, parameters are initialized
#' by uniformly segmenting the data into K segments, and for the next passes,
#' parameters are initialized by randomly segmenting the data into K contiguous
#'  segments.
#' @param max_iter The maximum number of iterations for the EM algorithm.
#' @param threshold A numeric value specifying the threshold for the relative
#'  difference of log-likelihood between two steps  of the EM as stopping
#'  criteria.
#' @param verbose A logical value indicating whether values of the
#' log-likelihood should be printed during EM iterations.
#' @param verbose_IRLS A logical value indicating whether values of the
#' criterion optimized by IRLS should be printed at each step of the EM
#' algorithm.
#' @return EM returns an object of class [ModelRHLP][ModelRHLP].
#' @seealso [ModelRHLP], [ParamRHLP], [StatRHLP]
#' @export
emRHLP <- function(X, Y, K, p, q = 1, variance_type = 2, n_tries = 1, max_iter = 1500, threshold = 1e-6, verbose = FALSE, verbose_IRLS = FALSE) {

  fData <- FData(X, Y)

  top <- 0
  try_EM <- 0
  best_loglik <- -Inf
  cpu_time_all <- c()

  while (try_EM < n_tries) {
    try_EM <- try_EM + 1
    message("EM try nr ", try_EM)
    time <- Sys.time()

    # Initialization
    param <- ParamRHLP$new(fData = fData, K = K, p = p, q = q, variance_type = variance_type)
    param$initParam(try_EM)
    iter <- 0
    converge <- FALSE
    prev_loglik <- -Inf

    stat <- StatRHLP(paramRHLP = param)

    while (!converge && (iter <= max_iter)) {
      stat$EStep(param)

      reg_irls <- param$MStep(stat, verbose_IRLS)
      stat$computeLikelihood(reg_irls)

      iter <- iter + 1
      if (verbose) {
        message("EM     : Iteration : ", iter, "  log-likelihood : "  , stat$loglik)
      }
      if (prev_loglik - stat$loglik > 1e-5) {
        message("!!!!! EM log-likelihood is decreasing from ", prev_loglik, "to ", stat$loglik)
        top <- top + 1
        if (top > 20)
          break
      }

      # Test of convergence
      converge <- abs((stat$loglik - prev_loglik) / prev_loglik) <= threshold
      if (is.na(converge)) {
        converge <- FALSE
      } # Basically for the first iteration when prev_loglik is Inf

      prev_loglik <- stat$loglik
      stat$stored_loglik[iter] <- stat$loglik
    } # End of the EM loop

    cpu_time_all[try_EM] <- Sys.time() - time

    if (stat$loglik > best_loglik) {
      statSolution <- stat$copy()
      paramSolution <- param$copy()
      if (param$K == 1) {
        statSolution$tau_ik <- matrix(stat$tau_ik, nrow = param$fData$m, ncol = 1)
        statSolution$pi_ik <- matrix(stat$pi_ik, nrow = param$fData$m, ncol = 1)
      }
      else{
        statSolution$tau_ik <- stat$tau_ik[1:param$fData$m, ]
        statSolution$pi_ik <- stat$pi_ik[1:param$fData$m, ]
      }

      best_loglik <- stat$loglik
    }
    if (n_tries > 1) {
      message("max value: ", stat$loglik)
    }
  }

  # Computation of c_ig the hard partition of the curves and klas
  statSolution$MAP()

  if (n_tries > 1) {
    message("max value: ", statSolution$loglik)
  }

  # Finish the computation of statistics
  statSolution$computeStats(paramSolution, cpu_time_all)

  return(ModelRHLP$new(paramRHLP = paramSolution, statRHLP = statSolution))
}
