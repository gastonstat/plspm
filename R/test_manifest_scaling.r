#' @title Test scaling of selected manifest variables
#' 
#' @details
#' Internal function. \code{test_manifest_scaling} is called by \code{plspm}.
#' Basically checks that unordered factors have nominal scaling.
#'
#' @param MV matrix or data frame with selected manifest variables
#' @param scaling list with specified scaling for each manifest variable
#' @return gives an error when factors in MV have bad scaling
#' @keywords internal
#' @template internals
#' @export
test_manifest_scaling <- function(MV, scaling)
{
  # to be used when MV is a data.frame containing factors
  if (is.data.frame(MV)) {
    mvs_class = sapply(MV, class)
    mvs_as_factors <- mvs_class == "factor"
    # if there are MVs as factors, check right scaling
    if (sum(mvs_as_factors) > 0) {
      factors_scaling = unlist(scaling)[mvs_as_factors]
      
      # factors can't be numeric or raw
      if (any(factors_scaling %in% c("num", "raw")))
        stop("\n'Data' contains factors that can't have metric scaling")
      
      # unordered factors must be nominal
      if (sum(factors_scaling == "ord") == 1) {
        if (!is.ordered(MV[,mvs_as_factors]))
          stop("\nUnordered factors in 'Data' can't have ordinal scaling")
      } 
      if (sum(factors_scaling == "ord") > 1)  {
        unordered = !apply(MV[,mvs_as_factors], 2, is.ordered)
        if (sum(unordered) > 0)
          stop("\nUnordered factors in 'Data' can't have ordinal scaling")
      }    
    }    
  }
  TRUE
}
