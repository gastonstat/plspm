#' @title Test outer weights convergence within specified maxiter
#' 
#' @details
#' Internal function. \code{test_weights_output} is called by \code{plspm}.
#' Basically checks that weights are not NULL.
#'
#' @param wgs output from 'get_weights()'
#' @param specs list of specifications
#' @return gives an error when wgs is NULL (i.e. non-convergence)
#' @keywords internal
#' @template internals
#' @export
test_null_weights <- function(wgs, specs) 
{
  if (is.null(wgs)) {
    print(paste("Iterative process is non-convergent with 'maxiter'=", 
                specs$maxiter, " and 'tol'=", specs$tol, sep=""))
    message("Algorithm stops") 
    stop("")
  }
  TRUE
}
