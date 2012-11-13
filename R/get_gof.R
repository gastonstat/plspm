#' @title Goodness-of-fit for \code{plspm}
#'
#' @description
#' Internal function. \code{get_gof} is called by \code{plspm}
#'
#' @param comu communalities
#' @param R2 R-squared coefficient
#' @param blocks list of variables in each block
#' @param IDM Inner Design Matrix
#' @keywords internal
#' @export
get_gof <-
function(comu, R2, blocks, IDM)
{
  lvs = nrow(IDM)
  blocklist = as.list(1:lvs)
  for (j in 1:lvs) blocklist[[j]] = rep(j, blocks[j])
  blocklist = unlist(blocklist)
  endo = rowSums(IDM)
  endo[endo != 0] = 1  
  n.end = sum(endo)
  # average of communalities
  R2.aux <- R2[endo == 1]
  comu.aux <- n.comu <- 0    
  for (j in 1:lvs)
  {
    if (length(which(blocklist==j)) > 1)
    {
      comu.aux = comu.aux + sum(comu[which(blocklist==j)])
      n.comu = n.comu + length(which(blocklist==j))
    }
  }
  gof = sqrt((comu.aux/n.comu) * mean(R2.aux))
  return(gof)
}
