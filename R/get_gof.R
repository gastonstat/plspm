#' @title Goodness-of-fit for \code{plspm}
#'
#' @details
#' Internal function. \code{get_gof} is called by \code{plspm}
#'
#' @param comu list of communalities
#' @param R2 vector of R-squared coefficients
#' @param blocks list of variables in each block
#' @param path_matrix Inner Design Matrix
#' @keywords internal
#' @template internals
#' @export get_gof
get_gof <- function(comu, R2, blocks, path_matrix)
{
  lvs = nrow(path_matrix)
  blocklist = indexify(blocks)  
  endo = rowSums(path_matrix)
  endo[endo != 0] = 1  
  
  # average of communalities
  R2_aux <- R2[endo == 1]
  comu_aux <- n_comu <- NULL
  for (j in 1:lvs)
  {
    # non mono factorial blocks only
    if (sum(blocklist==j) > 1)
    {
      comu_aux = c(comu_aux, mean(comu[blocklist==j]))
      n_comu = c(n_comu, sum(blocklist==j))
    }
  }
  mean_communality = sum(comu_aux * n_comu) / sum(n_comu)
  gof = sqrt(mean_communality * mean(R2_aux))
  # output
  gof
}
