#' @title Non-Metric Nominal Scale
#' 
#' @details
#' Internal function. \code{get_nom_scale} is called by \code{plspm}.
#'
#' @note
#' This function replaces the elements of x by the the means of y conditioned 
#' to the levels of x
#' 
#' @param y vector of values
#' @param x vector of the natural number series 1, 2, ..., p obtained 
#' by means of the function myRank
#' @param Xdummy dummy matrix corresponding to x
#' @return scaled matrix
#' @keywords internal
#' @template internals
#' @export
get_nom_scale <- function(y, x, Xdummy) 
{
  n = length(x)
  p = max(x, na.rm = TRUE)
  
  # vector of the means of y conditioned to the levels of x
  quant = tapply(y, x, mean, na.rm=TRUE)
  x_quant = Xdummy %*% quant
  x_quant
  
  # ===========  just in the case you need of them ============
  # ===========  eta2 = correlation rato  ============
  #eta2 <- var(x_quant)/var(y)
  #eta2<-var(x_quant,na.rm = T)/var(y, na.rm = T)	
  #list(Xdummy=Xdummy,x_quant=x_quant, eta2=eta2, quant=quant)
}
