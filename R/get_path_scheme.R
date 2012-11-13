#' @title Calculate inner weighting path scheme
#' 
#' @description
#' Internal function. \code{get_path_scheme} is called by \code{get_weights}
#' 
#' @param IDM Inner Design Matrix
#' @param Y Matrix of latent variables
#' @keywords internal
#' @export
get_path_scheme <-
  function(IDM, Y)
  {
    lvs <- nrow(IDM)
    E <- IDM
    for (k in 1:lvs) 
    {
      if (length(which(IDM[k,]==1)) > 0)
        E[which(IDM[k,]==1),k] <- lm(Y[,k]~Y[,which(IDM[k,]==1)]-1)$coef
      if (length(which(IDM[,k]==1)) > 0)
        E[which(IDM[,k]==1),k] <- cor(Y[,k], Y[,which(IDM[,k]==1)])
    }                 
    return(E)
  }
