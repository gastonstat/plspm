#' @title Calculate Dillon-Goldstein's rho
#' 
#' @details
#' Internal function. \code{get_rho} is called by \code{get_unidim} and 
#' \code{rho}.
#'
#' @param Block matrix (one block) of standardized manifest variables 
#' @param score first principal component of \code{Block} 
#' @return Dillon-Goldstein's rho
#' @keywords internal
#' @template internals
#' @export
get_rho <- function(Block, score)
{
  # dillon-goldstein rho
  p = ncol(Block)
  rho_numer = colSums(cor(Block, score))^2
  rho_denom = rho_numer + (p - colSums(cor(Block, score)^2) )
  Rho = rho_numer / rho_denom
  
  # output
  Rho
}
