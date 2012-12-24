#' @title Rescale Latent Variable Scores
#' 
#' @description 
#' Rescale standardized latent variable scores to original scale of manifest variables
#' 
#' @details
#' \code{rescale} requires all outer weights to be positive
#' 
#' @param pls object of class \code{"plspm"}
#' @param Y Optional dataset (matrix or data frame) used when argument
#' \code{dataset=NULL} inside \code{pls}.
#' @return A matrix with the rescaled latent variable scores
#' @author Gaston Sanchez
#' @seealso \code{\link{plspm}}
#' @export
#' @examples
#' 
#' \dontrun{
#'  ## example with customer satisfaction analysis
#'  
#'  # load data satisfaction
#'  data(satisfaction)
#'  
#'  # define inner model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0) 
#'  LOY = c(1,0,0,0,1,0)
#'  sat.inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'  
#'  # define outer model list
#'  sat.outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'  
#'  # define vector of reflective modes
#'  sat.mod = rep("A", 6)
#'  
#'  # apply plspm
#'  my_pls = plspm(satisfaction, sat.inner, sat.outer, sat.mod, scheme="factor", 
#'               scaled=FALSE)
#'               
#'  # rescaling standardized scores of latent variables
#'  scores = rescale(my_pls)
#'  
#'  # compare standardized LVs against rescaled LVs
#'  summary(my_pls$latents)
#'  summary(scores)
#'  }
#'
rescale <-
function(pls, Y = NULL)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (!inherits(pls, "plspm"))
    stop("\nSorry, an object of class 'plspm' was expected")
  # if Y available
  if (!is.null(Y))
  {
    if (is.null(pls$data))
    {
      if (!is.matrix(Y) && !is.data.frame(Y))
        stop("\n'Y' must be a numeric matrix or data frame")
      if (nrow(Y) != nrow(pls$latents))
        stop("\n'pls' and 'Y' are incompatible. Different number of rows")
    }
  } else { 
    # if no Y
    if (is.null(pls$data)) 
      stop("\n'Y' is missing; No dataset available in 'pls'")
  }
  # check positive outer weights
  wgs = pls$out.weights
  if (any(wgs < 0))
    stop("\nSorry, all outer weights must be positive")
  
  # =======================================================
  # prepare ingredients
  # =======================================================
  IDM = pls$model$IDM
  blocks = pls$model$blocks   
  modes = pls$model$modes
  lvs = nrow(IDM)
  mvs = sum(blocks)
  LVS = pls$latents
  outer = pls$model$outer
  
  # create block list
  blocklist = as.list(1:lvs)
  for (j in 1:lvs) 
    blocklist[[j]] = rep(j, blocks[j])
  blocklist = unlist(blocklist)

  # outer design matrix 'ODM' with normalized weights
  ODM = matrix(0, mvs, lvs)
  for (j in 1:lvs) 
    ODM[which(blocklist==j),j] = wgs[blocklist == j] / sum(wgs[blocklist == j])

  # calculating rescaled scores
  if (!is.null(pls$data)) 
  {
    # get rescaled scores
    Scores = pls$data %*% ODM
  } else {         
    # building data matrix 'DM'
    DM = matrix(NA, nrow(pls$latents), sum(blocks))
    for (k in 1:lvs)
      DM[,which(blocklist==k)] <- as.matrix(Y[,outer[[k]]])
    # get rescaled scores
    Scores = DM %*% ODM
  }

  # result
  dimnames(Scores) = list(rownames(pls$latents), rownames(IDM))
  Scores
}
