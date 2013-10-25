#' @title Non-Metric Numerical Scale
#' 
#' @details
#' Internal function. \code{get_num_scale} is called by \code{plspm}.
#'
#' @note
#' scales a matrix X in such a way that mean(X[,j])=0 and varpop(X[,j])=1
#' this means that sum(X[,j]^2) = n
#' if MD, sum(X[,j]^2, na.rm=T) = number of available elements
#' 
#' @param X a matrix to be scaled
#' @return scaled matrix
#' @keywords internal
#' @template internals
#' @export
get_num_scale <- function(X) {
  X = as.matrix(X)
  X_scaled = matrix(0, nrow(X), ncol(X))
  for (j in 1:ncol(X) ) {
    correction <- (sqrt(length(na.omit(X[,j]))/(length(na.omit(X[,j]))-1)))
    X_scaled[,j] <- scale(X[,j]) * correction
  }
  #rownames(X_scaled) = rownames(X) 
  #colnames(X_scaled) = colnames(X)
  X_scaled
}
