source("R/dataset.R")

MixModel <- setRefClass(
  "MixModel",
  contains = "MyData",
  # Define the fields
  fields = list(K="numeric", # nombre de regimes
                p="numeric", # dimension de beta (ordre de reg polynomiale)
                q="numeric" # dimension de w (ordre de reg logistique)
                ),
  methods = list(

  )
)

MixModel<-function(mixData, K,p,q){
  new("MixModel",y=mixData$y, x=mixData$x, m=mixData$m, K=K, p=p, q=q)
}


