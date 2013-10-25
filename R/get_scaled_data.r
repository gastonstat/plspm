#' @title Scaling data outside \code{plspm}
#' 
#' @details
#' Internal function. \code{get_scaled_data} is called by \code{plspm.groups}
#' and others.
#' 
#' @note 
#' In non-metric case, all mvs are starndardized, and those ordinal or nominal
#' are rankified
#' 
#' @param DM matrix or data frame from with manifest variables
#' @param scaled whether to scale latent variables
#' @return matrix or data frame of (un)scaled MVs 
#' @keywords internal
#' @template internals
#' @export
get_data_scaled <- function(DM, scaled) {
  if (scaled) {
    sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
    X <- scale(DM, scale=sd.X)
  } else {
    X <- scale(DM, scale=FALSE)
  }
  # output
  X
}
