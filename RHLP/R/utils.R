normalize <- function(A, dim){
# NORMALIZE makes the entries of a (multidimensional <= 2) array sum to 1.

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
