#' @title Rank of a non-metric variable
#' 
#' @details
#' Internal function. \code{get_rank} is called by \code{get_treated_data}
#' and \code{get_weights}.
#'
#' @param X an ordinal or nominal manifest variable to be ranked
#' @return vector with the corresponding rank for each category
#' @keywords internal
#' @template internals
#' @export
get_rank <- function(X) 
{
  X_ranked = rep(0, length(X))
  uniq = unique(X, ties='min')
  uniq_ranked = rank(uniq)
  for (k in 1:length(uniq)) {
    X_ranked[which(X == uniq[k])] <- uniq_ranked[k]  
  }

  X_ranked
}
