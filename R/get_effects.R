#' @title Path coefficient effects for \code{plspm}
#' 
#' @details
#' Internal function. \code{get_effects} is called by \code{plspm}.
#' 
#' @param Path matrix of path coefficients
#' @return data frame with relationships, direct, indirect, and total effects
#' @keywords internal
#' @template internals
#' @export
get_effects <- function(Path)
{
  # how many latent variables and their names
  lvs = nrow(Path)
  lvs_names = rownames(Path)
  
  # list for storing effects
  path_effects = as.list(1:(lvs-1))
  path_effects[[1]] = Path
  
  # when only 2 lvs
  if (lvs == 2) {
    indirect_paths = matrix(c(0,0,0,0), 2, 2)
    total_paths = Path
  } else {
    # when more than 2 lvs
    for (k in 2:(lvs-1)) {
      path_effects[[k]] = path_effects[[k-1]] %*% Path        
    }
    indirect_paths = matrix(0, lvs, lvs)
    for (k in 2:length(path_effects)) {
      indirect_paths = indirect_paths + path_effects[[k]]        
    }
    total_paths = Path + indirect_paths
  }
  
  # initialize
  efs_labels <- direct <- indirect <- total <- NULL
  for (j in 1:lvs) {
    for (i in j:lvs) {
      if (i != j) {
        efs_labels = c(efs_labels, paste(lvs_names[j], "->", lvs_names[i]))
        direct = c(direct, Path[i,j])
        indirect = c(indirect, indirect_paths[i,j])
        total = c(total, total_paths[i,j])
      }
    }
  }
  
  # output
  data.frame(relationships = efs_labels, 
             direct = direct, 
             indirect = indirect, 
             total = total)
}
