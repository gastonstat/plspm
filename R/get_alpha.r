#' @title Calculate Cronbach's alpha
#' 
#' @details
#' Internal function. \code{get_alpha} is called by \code{get_unidim} and 
#' \code{alpha}.
#'
#' @param Block matrix (one block) of standardized manifest variables 
#' @return Cronbach's alpha
#' @keywords internal
#' @template internals
#' @export
get_alpha <- function(Block)
{
  # how many variables
  p = ncol(Block)
  correction = sqrt((nrow(Block)-1) / nrow(Block)) 
  
  # cronbach's alpha
  alpha_denom = var(rowSums(Block)) * correction^2
  alpha_numer = 2 * sum(cor(Block)[lower.tri(cor(Block))])
  alpha = (alpha_numer / alpha_denom) * (p / (p - 1))
  Alpha <- ifelse(alpha < 0, 0, alpha)
  
  # output
  Alpha
}
