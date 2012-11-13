#' @title Path coefficient effects for \code{plspm}
#' 
#' @description
#' Internal function. \code{get_effects} is called by \code{plspm}.
#' 
#' @param Path Matrix of path coefficients
#' @return data frame with direct, indirect, and total effects
#' @keywords internal
#' @export
get_effects <-
function(Path)
{
  # how many latent variables and their names
  lvs = nrow(Path)
  lvs.names = rownames(Path)
  # list to store effects
  path.efects = as.list(1:(lvs-1))
  path.efects[[1]] = Path
  # when only 2 lvs
  if (lvs == 2) {
    ind.paths = matrix(c(0,0,0,0), 2, 2)
    total.paths = Path
  } else {
    # when more than 2 lvs
    for (k in 2:(lvs-1))
      path.efects[[k]] = path.efects[[k-1]] %*% Path
    # indirect effects
    ind.paths = matrix(0, lvs, lvs)
    for (k in 2:length(path.efects))
      ind.paths = ind.paths + path.efects[[k]]
    # total effects
    total.paths = Path + ind.paths
  }
  # initialize
  efs.labs <- dir.efs <- ind.efs <- tot.efs <- NULL
  for (j in 1:lvs) 
  {
    for (i in j:lvs)
    {
      if (i != j) 
      {
        # labels
        efs.labs = c(efs.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
        # direct effects
        dir.efs = c(dir.efs, Path[i,j])
        # indirect effects
        ind.efs = c(ind.efs, ind.paths[i,j])
        # total effects
        tot.efs = c(tot.efs, total.paths[i,j])
      }
    }
  }
  # results
  Effects = data.frame(relationships = efs.labs, 
                       dir.effects = dir.efs, 
                       ind.effects = ind.efs, 
                       tot.effects = tot.efs)
  return(Effects)
}
