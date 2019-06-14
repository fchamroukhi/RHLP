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
