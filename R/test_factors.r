#' @title Test presence of factors
#' 
#' @details
#' Internal function. \code{test_factors} is called by 
#' \code{get_treated_data}.
#' 
#' @note Test presence of factors in a data.frame
#' 
#' @param DF data frame that might contain factors
#' @return whether DF contains factors or not
#' @keywords internal
#' @template internals
#' @export
test_factors <- function(DF)
{
  contains_factors = FALSE
  # to be used when DF is a data.frame containing factors
  if (is.data.frame(DF)) {
    mvs_class = sapply(DF, class)
    mvs_as_factors <- mvs_class == "factor"
    # tell me if there are MVs as factors
    if (sum(mvs_as_factors) > 0)
      contains_factors = TRUE
  }
  contains_factors
}


#' @title Transform factors in MV into numeric
#' 
#' @details
#' Internal function. \code{get_numerics} is called by 
#' \code{get_treated_data}.
#' 
#' @note This function is used when there are categorical variables 
#' defined as 'factors'. The idea is to convert those factors into numeric
#' vectors while keeping the levels (ie categories) in a separate list
#' 
#' @param MV data frame with manifest variables as factors
#' @return list with transformed MV and categories of factors
#' @keywords internal
#' @template internals
#' @export
get_numerics <- function(MV)
{  
  mvs_class = sapply(MV, class)
  mvs_as_factors <- mvs_class == "factor"
  categorical = which(mvs_as_factors)
  categories = vector(mode="list", ncol(MV))
  
  # only keep levels of categorical mvs in 'factor' format
  for (j in seq_along(categorical)) {
    # keep original levels
    categories[[categorical[j]]] = levels(MV[,categorical[j]])
    # convert to numeric
    MV[,categorical[j]] = as.numeric(MV[,categorical[j]])
  }    
  
  # output
  list(MV=MV, categories=categories)
}
