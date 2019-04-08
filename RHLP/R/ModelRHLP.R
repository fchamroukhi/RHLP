source("R/MetaParamRHLP.R")
source("R/ParamRHLP.R")
source("R/StatRHLP.R")

ModelRHLP <- setRefClass(
  "ModelRHLP",
  contains = c("MetaParamRHLP", "ParamRHLP")

)

ModelRHLP<-function(K, p, q, variance_type){
  metap <- MetaParamRHLP(K, p, q, variance_type)
  paramp <- ParamRHLP(metap)
  new("ModelRHLP", K = metap$K, p = metap$p, q = metap$q, variance_type = metap$variance_type, W = paramp$W, beta = paramp$beta, sigma = paramp$sigma)
}
