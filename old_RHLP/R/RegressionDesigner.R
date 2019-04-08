RegressionDesigner <- setRefClass(
  "RegressionDesigner",
  fields = list(
    XBeta = "matrix",
    Xw = "matrix"
  ),

  methods = list(
    ##pour 1 courbe
    designmatrix_FRHLP = function(x,p,q=NULL){
      if (ncol(x) == 1){
        x<-t(x)
      }
      order_max <- p
      if (!is.null(q)){
        order_max <- max(p,q)
      }

      X <- matrix(NA, length(x), order_max+1)
      for (i in 0:(order_max)){
        X[,i+1] <- x^i # phi2w = [1 t t.^2 t.^3 t.^p;......;...]
      }

      XBeta <<- X[,1:(p+1)]; # Matrice de regresseurs pour Beta
      if (!is.null(q)){
        Xw <<- X[,1:(q+1)]; # matrice de regresseurs pour w
      }
    },

    setPhi1 = function(x,p,q){
      ##pour 1 courbe
      designmatrix_FRHLP(x, p, q)
    }
  )
)
