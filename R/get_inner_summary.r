#' @title Inner summary assessment
#' 
#' @details
#' Internal function. \code{get_inner_summary} is called by \code{plspm}.
#'
#' @param path_matrix matrix of path connections
#' @param blocks list indicating blocks of manifest variables
#' @param modes vector of modes
#' @param communality vector with communality values
#' @param redundancy vector with redundancy values
#' @param R2 vector with R2 values
#' @return A data frame with the following columns:
#' @return \item{Type}{Exogenous or Endogenous}
#' @return \item{R2}{R2 coefficient}
#' @return \item{Mean_Communality}{average communality}
#' @return \item{Mean_Redundancy}{average redudancy}
#' @return \item{AVE}{Average Variance Extracted}
#' @keywords internal
#' @template internals
#' @export
get_inner_summary <- 
function(path_matrix, blocks, modes, communality, redundancy, R2)
{
  blocklist = indexify(blocks)  
  exo_endo = rep("Exogenous", nrow(path_matrix))
  exo_endo[rowSums(path_matrix) != 0] = "Endogenous"
  avg_comu = rep(0, nrow(path_matrix))
  avg_redu = rep(0, nrow(path_matrix))
  AVE = rep(0, nrow(path_matrix))
  
  for (k in 1:nrow(path_matrix))
  {
    avg_comu[k] = mean(communality[blocklist == k])
    avg_redu[k] = mean(redundancy[blocklist == k])
    if (modes[k] == "A")
    {
      ave_num = sum(communality[blocklist == k])
      ave_denom = ave_num + sum(1 - (communality[blocklist == k]))
      AVE[k] = ave_num / ave_denom
    }
  }
  
  # output
  data.frame(Type = exo_endo, 
             R2 = R2, 
             Block_Communality = avg_comu, 
             Mean_Redundancy = avg_redu, 
             AVE = AVE,
             row.names = rownames(path_matrix))
}
