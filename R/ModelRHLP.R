#' A Reference Class which represents a RHLP model.
#'
#' ModelRHLP represents a RHLP model. It inherits from [FData][FData] class. It
#' contains all meta-paramters such as `K`, the number of segments, `p` the
#' order of the polynomial regression\dots
#'
#' @usage NULL
#' @field K The number of regimes (mixture components).
#' @field p The order of the polynomial regression.
#' @field q The dimension of the logistic regression. For the purpose of
#' segmentation, it must be set to 1.
#' @field variance_type Numeric indicating if the model is homoskedastic
#' (`variance_type` = 1) or heteroskedastic (`variance_type` = 2).
#' @seealso [FData]
#' @export
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
    nu <- (p + q + 3) * K - (q + 1) - (K - 1)
  } else{
    nu <- (p + q + 3) * K - (q + 1)
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
