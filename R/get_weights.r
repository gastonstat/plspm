#' @title Outer Weights
#' 
#' @details
#' Internal function. \code{get_weights} is called by \code{plspm}
#' 
#' @note
#' Calculate outer weights (under Lohmoller's algorithm)
#' 
#' @param X scaled data
#' @param path_matrix matrix with path connections
#' @param blocks list with variables in each block
#' @param specs list with algorithm specifications
#' @return list of outer weights, ODM, iter
#' @export
#' @template internals
#' @keywords internal
get_weights <- function(X, path_matrix, blocks, specs)
{
  lvs = nrow(path_matrix)
  mvs = ncol(X)
  sdv = sqrt((nrow(X)-1) / nrow(X))   # std.dev factor correction
  blockinds = indexify(blocks)

  # outer design matrix 'ODM' and matrix of outer weights 'W'
  ODM = list_to_dummy(blocks)
  W = ODM %*% diag(1/(apply(X %*% ODM, 2, sd)*sdv), lvs, lvs)
  w_old = rowSums(W)    
  iter = 1
  
  repeat 
  {
    # external estimation of LVs 'Y'
    Y = X %*% W
    Y = scale(Y) * sdv
    # matrix of inner weights 'e' 
    E <- switch(specs$scheme, 
                "centroid" = sign(cor(Y) * (path_matrix + t(path_matrix))),
                "factorial" = cor(Y) * (path_matrix + t(path_matrix)),
                "path" = get_path_scheme(path_matrix, Y))
    # internal estimation of LVs 'Z'
    Z = Y %*% E  
#    Z = Z %*% diag(1/(apply(Z,2,sd)*sdv), lvs, lvs)  # scaling Z
    # computing outer weights 'w'
    for (j in 1:lvs)
    {
      if (specs$modes[j] == "A")
        W[blockinds==j,j] <- (1/nrow(X)) * Z[,j] %*% X[,blockinds==j] 
      if (specs$modes[j] == "B")
        W[blockinds==j,j] <- solve.qr(qr(X[,blockinds==j]), Z[,j])
    }
    w_new = rowSums(W)                
    w_dif = sum((abs(w_old) - abs(w_new))^2) 
    if (w_dif < specs$tol || iter == specs$maxiter) break
    w_old = w_new
    iter = iter + 1
  } # end repeat       
    
  # preparing results
  if (iter == specs$maxiter) {
    results = NULL
  } else {
    W = W %*% diag(1/(apply(X %*% W, 2, sd)*sdv), lvs, lvs)
    w_new = rowSums(W)                
    names(w_new) = colnames(X)
    dimnames(W) = list(colnames(X), rownames(path_matrix))    
    results = list(w = w_new, W = W, ODM = ODM, iter = iter)
  }
  # output
  results
}
