#' @title Calculate AVE (Average Variance Extracted)
#' 
#' @details
#' Internal function. \code{get_ave} is called by \code{plspm}.
#' 
#' @param communality list of communalities for each block
#' @return vector of average variance extracted for each block
#' @keywords internal
#' @template internals
#' @export
get_ave <- function(communality)
{
  # average variance extracted
  sapply(communality, function(x) sum(x) / (sum(x) + sum(1-x)))
}
