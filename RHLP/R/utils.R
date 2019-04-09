source("R/FData.R")


ones <- function(n,d,g=1){
  if (g==1){
    return(matrix(1,n,d))
  }
  else{
    return(array(1,dim=c(n,d,g)))
  }
}

zeros <- function(n,d,g=1){
  if (g==1){
    return(matrix(0,n,d))
  }
  else{
    return(array(0,dim=c(n,d,g)))
  }
}

rand <- function(n,d,g=1){
  if (g==1){
    return(matrix(runif(n*d), n,d))
  }
  else{
    return(array(runif(n*d),dim=c(n,d,g)))
  }
}

repmat <- function(M, n, d){
  return(kronecker(matrix(1, n, d), M))
}

drnorm <- function(n, d, mean, sd){
  A <- matrix(nrow = n, ncol = d)
  for (i in 1:d){
    A[,i] <- rnorm(n,mean,sd)
  }
  return(A)
}


lognormalize <- function(M){
  if (!is.matrix(M)){
    M <- matrix(M)
  }
  n <- nrow(M)
  d <- ncol(M)
  a <- apply(M,1,max)
  return(M-repmat(a + log(rowSums(exp(M - repmat(a,1,d)))), 1, d))
}

normalize <- function(A, dim){
  # Normalize makes the entries of a (multidimensional <= 2) array sum to 1.
  # Input
  # A = Array to be normalized
  # dim = dimension is specified to normalize.
  # Output
  # M = Array after normalize.
  # z is the normalize constant
  # Note:
  # If dim is specified, we normalize the specified dimension only,
  # Otherwise we normalize the whole array.
  # Dim = 1 normalize each column
  # Dim = 2 normalize each row

  if (nargs() < 2 ){
    z = sum(A)
    # Set any zeros to one before dividing
    # This is valid, since c = 0 ==> all i.A[i] = 0 ==> the anser should be 0/1 = 0.
    s = z + (z==0)
    M = A / s
  } else if (dim == 1){ # normalize each column
    z = colSums(A)
    s = z + (z==0)
    M = A / matrix(s, nrow=dim(A)[1], ncol=length(s),byrow=TRUE)
  } else{
    z = rowSums(A)
    s = z + (z==0)
    M = A / matrix(s, ncol=dim(A)[2], nrow=length(s),byrow = FALSE)
  }
  output = list("M" = M, "z" = z)
  return(output)
}















generateRandomDataSet = function(){
  nn<-c(10,10,10)
  n1<-nn[1]
  n2<-nn[2]
  n3<-nn[3]
  mean_flou <- as.matrix(read.table("data/mean_1_flou.txt",
                                    header = FALSE))
  y1 <- ones(n1,1)%*%t(mean_flou) + t(drnorm(length(mean_flou),n1,5,1))+1
  a <- drnorm(80,n2,7,1)
  b <- drnorm(130,n2,5,1)
  c <- drnorm(140,n2,4,1)
  y3 <- t(rbind(a, b, c))

  a <- drnorm(120,n3,5,1)
  b <- drnorm(70,n3,7,1)
  c <- drnorm(160,n3,5,1)
  y2 <- t(rbind(a, b, c))
  Y <- rbind(y1, y2, y3)

  m <- nrow(Y)
  X <- t(matrix(seq(0, 1, length.out = m)))

  fData <- FData$new()
  fData$setData(X, Y)
  return(fData)
}

#d <- generateRandomDataSet()
