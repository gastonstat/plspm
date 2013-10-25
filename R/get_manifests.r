#' @title Building data matrix with manifest variables
#' 
#' @details
#' Internal function. \code{get_manifests} is called by \code{plspm}.
#' 
#' @param Data matrix or data frame from which manifest variables are extracted
#' @param blocks list (wtih numeric elements) indicating the set of manifest 
#' variables that form each block
#' @return matrix or data frame of selected manifest variables
#' @keywords internal
#' @template internals
#' @export
get_manifests <- function(Data, blocks)
{
  # building data matrix 'MV'
  MV = Data[,unlist(blocks)]
  
  # add row and column names
  mvs_names = colnames(Data)[unlist(blocks)]
  dimnames(MV) = list(rownames(Data), mvs_names)

  # output
  MV
}
