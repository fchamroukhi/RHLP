ModelRHLP <- setRefClass(
  "ModelRHLP",
  fields = list(K="numeric", # nombre de regimes
                p="numeric", # dimension de beta (ordre de reg polynomiale)
                q="numeric", # dimension de w (ordre de reg logistique)
                variance_type="numeric",
                W = "matrix",
                beta = "matrix",
                sigma = "matrix"
  )

)

ModelRHLP<-function(K, p, q, variance_type){
  W <- matrix(0,q+1, K-1)
  beta <- matrix(NA, p+1, K)
  if (variance_type == variance_types$homoskedastic){
    # todo: verify the dimensions
    sigma <- matrix(1)
  }
  else{
    sigma <- matrix(NA, K)
  }

  new("ModelRHLP", K = K, p = p, q = q, variance_type = variance_type, W = W, beta = beta, sigma = sigma)
}
