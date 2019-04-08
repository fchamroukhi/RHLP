source("R/enums.R")
source("R/utils.R")

ParamRHLP <- setRefClass(
  "ParamRHLP",
  fields = list(
    W = "matrix",
    beta = "matrix",
    sigma = "matrix"
  )
)

ParamRHLP<-function(metaParamRHLP){
  W <- matrix(0,metaParamRHLP$q+1, metaParamRHLP$K-1)
  beta <- matrix(NA, metaParamRHLP$p+1, metaParamRHLP$K)
  if (metaParamRHLP$variance_type == variance_types$homoskedastic){
    # todo: verify the dimensions
    sigma <- matrix(NA, metaParamRHLP$K)
  }
  else{
    sigma <- matrix(NA, metaParamRHLP$K)
  }
  new("ParamRHLP", W = W, beta = beta, sigma = sigma)
}
