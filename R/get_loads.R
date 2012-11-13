#' @title Calculate loadings and communality for \code{plspm}
#' 
#' @description
#' Internal function. \code{get_loads} is called by \code{plspm}.
#' 
#' @param X Data matrix
#' @param Y.lvs Matrix of latent variables
#' @param blocks vector with number of manifest variabls per block
#' @return oadinsg and communalities
#' @keywords internal
#' @export
get_loads <-
  function(X, Y.lvs, blocks)
  {
    lvs = length(blocks)
    mvs = ncol(X)
    blocklist = as.list(1:lvs)
    for (j in 1:lvs) blocklist[[j]] = rep(j, blocks[j])
    blocklist = unlist(blocklist)
    loads <- comu <- rep(NA, mvs)    
    for (j in 1:lvs)
      loads[blocklist==j] = cor(X[,blocklist==j], Y.lvs[,j])
    comu = loads^2
    names(loads) <- names(comu) <- colnames(X)
    res.loads = list(loads, comu)
    return(res.loads)
  }
