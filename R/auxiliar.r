#' @title Presence of missing values 
#' 
#' @details
#' Internal function. \code{is_missing} is called by 
#' \code{get_weights_nonmetric}
#'
#' @param x a list
#' @return vector indicating whether each element in x contains NAs
#' @keywords internal
#' @template internals
#' @export
is_missing <- function(x) 
{
  if (is.matrix(x)) {
    num_miss = apply(x, 2, function(u) sum(is.na(u)))    
  } else {
    num_miss = sum(is.na(x))
  }
  as.logical(sum(num_miss))
}


#' @title Normalize a vector 
#' 
#' @details
#' Internal function. \code{normalize} is called internally
#'
#' @param x a numeric vector
#' @return normalized vector
#' @keywords internal
#' @template internals
#' @export
normalize <- function(x)
{
  if (!is.vector(x) || !is.numeric(x))
    stop("\nA numeric vector is required") 
  # output
  x / sqrt(sum(x * x))
}
