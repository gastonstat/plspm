#' @title Test Data Set Availibility
#' 
#' @details
#' Internal function. \code{test_dataset} checks whether a dataset is
#' available for plspm.groups, rescale, it.reb, rebus.pls, rebus.test
#' 
#' @param Dataset optional data set (with manifest variables) provided by user
#' @param pls_data Data matrix containing the manifest variables used in the
#' model. Only available when \code{dataset=TRUE} inside \code{plspm()}
#' @param num_obs number of rows in PLS-PM Scores
#' @return TRUE if dataset is available, otherwise FALSE
#' @keywords internal
#' @template internals
#' @export
test_dataset <- function(Dataset, pls_data, num_obs)
{
  if (!is.null(Dataset)) # if Y available
  {
    if (is.null(pls_data))
    {
      if (!is.matrix(Dataset) && !is.data.frame(Dataset))
        stop("Invalid object 'Y'. Must be a numeric matrix or data frame.")
      if (nrow(Dataset) != num_obs)
        stop("Argument 'pls' and 'Y' are incompatible. Different number of rows.")
    }
  } else { # if no Y
    if (is.null(pls_data)) 
      stop("Argument 'Y' is missing. No dataset available.")
  }
  # otherwise
  TRUE
}
