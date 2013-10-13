#' @title Get general parameters
#' 
#' @details
#' Internal function. \code{get_generals} is called by \code{plspm}.
#' This function gets the number and names of: observations,
#' manifest variables, and latent variables
#'
#' @param MV matrix of manifest variables
#' @param path_matrix matrix with path connections
#' @return list with number and names of observations, MVs and LVs
#' @keywords internal
#' @template internals
#' @export
get_generals <- function(MV, path_matrix)
{
  list(obs = nrow(MV),
       obs_names = rownames(MV),
       mvs = ncol(MV),
       mvs_names = colnames(MV),
       lvs = nrow(path_matrix),
       lvs_names = rownames(path_matrix))
}
