designmatrix = function(x, p, q = NULL, n=1) {
  order_max <- p
  if (!is.null(q)) {
    order_max <- max(p, q)
  }

  X <- matrix(NA, length(x), order_max + 1)
  for (i in 0:(order_max)) {
    X[, i + 1] <- x ^ i
  }

  XBeta <-X[, 1:(p + 1)]
  # design matrix for Beta (the polynomial regressions)
  if (!is.null(q)) {
    Xw <- X[, 1:(q + 1)]
    # design matrix for w (the logistic regression)
  }

  XBeta <- repmat(XBeta, n, 1)
  Xw <- repmat(Xw, n, 1)

  return(list(Xw, XBeta))
}
