FData <- setRefClass(
  "FData",
  fields = list(
    X = "numeric", # Covariates
    Y = "matrix", # Response
    m = "numeric",
    n = "numeric",
    vecY = "matrix"
  )
)

FData <- function(X, Y) {

  Y <- as.matrix(Y)

  n <- nrow(Y)
  m <- ncol(Y)

  vecY <- matrix(t(Y), ncol = 1)

  if (n == 1) {
    Y <- t(Y)
  }

  fData <- new("FData", X = X, Y = Y, m = m, n = n, vecY = vecY)
}
