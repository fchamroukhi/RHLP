source("R/dataset.R")

MetaParamRHLP <- setRefClass(
  "MetaParamRHLP",
  # Define the fields
  fields = list(K="numeric", # nombre de regimes
                p="numeric", # dimension de beta (ordre de reg polynomiale)
                q="numeric", # dimension de w (ordre de reg logistique)
                variance_type="numeric"
                )
)

MetaParamRHLP<-function(K, p, q, variance_type){
  new("MetaParamRHLP", K = K, p = p, q = q, variance_type = variance_type)
}


