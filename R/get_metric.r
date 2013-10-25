#' @title Type of metric based on scaling measurement
#' 
#' @details
#' Internal function. \code{get_metric} is called by \code{plspm}.
#' It is used to decide what function apply for the iterative process of the
#' PLS-PM algorithm.
#'
#' @note If scaling is NULL, apply get_weights(),
#' otherwise apply get_weights_nonmetric()
#' @param scaling list with measurement scale of each manifest variable
#' @return metric type of metric (metric data = TRUE / non-metric data = FALSE)
#' @keywords internal
#' @template internals
#' @export
get_metric <- function(scaling)
{
  # metric case
  if (is.null(scaling)) metric = TRUE else metric = FALSE
  
  # output
  metric
}
