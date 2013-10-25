#' @title Calculate inner weighting path scheme
#' 
#' @details
#' Internal function. \code{get_path_scheme} is called by \code{get_weights}
#' 
#' @param path_matrix matrix with paths
#' @param LV matrix of latent variables
#' @keywords internal
#' @template internals
#' @export
get_path_scheme <- function(path_matrix, LV)
{
  # matrix for inner weights
  E = path_matrix
  
  for (k in 1:ncol(path_matrix)) 
  {
    # followers
    follow <- path_matrix[k,] == 1
    if (sum(follow) > 0)
      E[follow,k] <- lm(LV[,k] ~ LV[,follow] - 1)$coef
    # predecesors
    predec <- path_matrix[,k] == 1
    if (sum(predec) > 0)
      E[predec,k] <- cor(LV[,k], LV[,predec])
  } 
  
  # output
  E
}
