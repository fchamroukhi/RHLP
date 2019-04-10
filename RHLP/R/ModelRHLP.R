source("R/FData.R")
source("R/enums.R")

ModelRHLP <- setRefClass(
  "ModelRHLP",
  contains = "FData",
  # Define the fields
  fields = list(K = "numeric", # number of regimes
                p = "numeric", # dimension of beta (order of polynomial regression)
                q = "numeric", # dimension of w (order of logistic regression)
                variance_type = "numeric",
                nu = "numeric" # degree of freedom
                )
)

ModelRHLP <- function(fData, K, p, q){

  if (variance_type == variance_types$homoskedastic){
    nu <<- (p + q + 3) * K - (q + 1) - (K - 1)
  }
  else{
    nu <<- (p + q + 3) * K - (q + 1)
  }

  new("ModelRHLP", Y = fData$Y, X = fData$X, m = fData$m, K = K, p = p, q = q, variance_type = variance_type, nu = nu)
}


