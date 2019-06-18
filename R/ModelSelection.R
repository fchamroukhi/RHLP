#' selectRHLP is used to select a RHLP model based on criteria such as BIC or AIC.
#'
#' @details selectRHLP is used to select the "best" RHLP model according to some
#' criteria such as BIC or AIC. This function runs every RHLP model by varying
#' the number of regimes `K` from `Kmin` to `Kmax` and the order of the
#' polynomial regression `p` from `pmin` to `pmax`.
#'
#' @param X Numeric vector of length \emph{m} representing the covariates.
#' @param Y Matrix of size \eqn{(n, m)} representing \emph{n} functions of `X`
#' observed at points \eqn{1,\dots,m}.
#' @param Kmin The minimum number of regimes (mixture components).
#' @param Kmax The maximum number of regimes (mixture components).
#' @param pmin The minimum order of the polynomial regression.
#' @param pmax The maximum order of the polynomial regression.
#' @param criterion The criterion used to select the "best" RHLP model.
#' @return selectRHLP returns an object of class [ModelRHLP][ModelRHLP]
#' representing the "best" RHLP model according to the selected `criterion`.
#' @seealso [ModelRHLP]
#' @export
selectRHLP <- function(X, Y, Kmin = 1, Kmax = 10, pmin = 0, pmax = 4, criterion = c("BIC", "AIC")) {

  criterion <- match.arg(criterion)

  vrhlp <- Vectorize(function(K, p, X1 = X, Y1 = Y) emRHLP(X = X1, Y = Y1, K, p),
                     vectorize.args = c("K", "p"))

  rhlp <- outer(Kmin:Kmax, pmin:pmax, vrhlp)

  if (criterion == "BIC") {
    results <- apply(rhlp, 1:2, function(x) x[[1]]$statRHLP$BIC)
  } else {
    results <- apply(rhlp, 1:2, function(x) x[[1]]$statRHLP$AIC)
  }
  rownames(results) <- sapply(Kmin:Kmax, function(x) paste0("(K = ", x, ")"))
  colnames(results) <- sapply(pmin:pmax, function(x) paste0("(p = ", x, ")"))


  selected <- rhlp[which(results == max(results), arr.ind = T)][[1]]

  cat(paste0("The RHLP model selected via the \"", criterion, "\" has K = ",
             selected$paramRHLP$K, " regimes \n and the order of the ",
             "polynomial regression is p = ", selected$paramRHLP$p, "."))
  cat("\n")
  cat(paste0("BIC = ", selected$statRHLP$BIC, "\n"))
  cat(paste0("AIC = ", selected$statRHLP$AIC, "\n"))

  return(selected)

}
