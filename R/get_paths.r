#' @title Calculate path coefficients for \code{plspm}
#' 
#' @details
#' Internal function. \code{get_paths} is called by \code{plspm}.
#' 
#' @param path_matrix path matrix
#' @param Y_lvs Matrix of latent variables
#' @param full logical to indicate all results from 'summary(lm())'
#' @return list with inner results, path coefs matrix, R2, and residuals
#' @keywords internal
#' @template internals
#' @export
get_paths <-  function(path_matrix, Y_lvs, full=TRUE)
{
  lvs_names = colnames(path_matrix)
  endogenous = as.logical(rowSums(path_matrix))
  num_endo = sum(endogenous)
  results = as.list(1:num_endo)
  Path = path_matrix
  residuals = as.list(1:num_endo)
  R2 = rep(0, nrow(path_matrix))
  
  for (aux in 1:num_endo) 
  {
    # index for endo LV
    k1 <- which(endogenous)[aux]
    # index for indep LVs
    k2 = which(path_matrix[k1,] == 1)
    
    path_lm = summary(lm(Y_lvs[,k1] ~ Y_lvs[,k2]))
    Path[k1,k2] = path_lm$coef[-1,1]
    residuals[[aux]] = path_lm$residuals  
    R2[k1] = path_lm$r.squared
    inn_val = c(path_lm$r.squared, path_lm$coef[,1])
    # ----- NEW results
    inn_labels = c("Intercept", names(k2))
    rownames(path_lm$coefficients) = inn_labels
    results[[aux]] <- path_lm$coefficients
    # ----- OLD results
    # inn_lab = c("R2", "Intercept", 
    # paste(rep("path_",length(k2)),names(k2),sep=""))
    # names(inn_val) = NULL
    # results[[aux]] <- data.frame(concept=inn_lab, value=round(inn_val, 4))
  }
  names(results) = lvs_names[endogenous]  
  names(R2) = lvs_names
  
  # output
  list(results, Path, R2, residuals)
}
