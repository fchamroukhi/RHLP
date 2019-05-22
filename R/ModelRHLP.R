ModelRHLP <- setRefClass(
  "ModelRHLP",
  contains = "FData",
  fields = list(
    K = "numeric", # Number of regimes
    p = "numeric", # Dimension of beta (order of polynomial regression)
    q = "numeric", # Dimension of w (order of logistic regression)
    variance_type = "numeric",
    nu = "numeric" # Degree of freedom
  )
)

ModelRHLP <- function(fData, K, p, q, variance_type) {
  if (variance_type == variance_types$homoskedastic) {
    nu <<- (p + q + 3) * K - (q + 1) - (K - 1)
  } else{
    nu <<- (p + q + 3) * K - (q + 1)
  }

  new(
    "ModelRHLP",
    Y = fData$Y,
    X = fData$X,
    m = fData$m,
    n = fData$n,
    K = K,
    p = p,
    q = q,
    variance_type = variance_type,
    nu = nu
  )
}
