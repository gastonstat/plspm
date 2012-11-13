#' @title Calculate outer weights for \code{plspm}
#' 
#' @description
#' Internal function. \code{get_weights} is called by \code{plspm}
#' 
#' @param X scaled data
#' @param IDM Inner Design Matrix
#' @param blocks list with variables in each block
#' @param modes vector measurement mode
#' @param scheme inner path scheme to use
#' @param tol tolerance threshold for convergen
#' @param iter maximum number of iterations
#' @export
#' @keywords internal
get_weights <-
function(X, IDM, blocks, modes, scheme, tol, iter)
{
  lvs = nrow(IDM)
  mvs = ncol(X)
  sdv = sqrt((nrow(X)-1) / nrow(X))   # std.dev factor correction
  blocklist <- as.list(1:lvs)
  for (j in 1:lvs) blocklist[[j]] = rep(j, blocks[j])
  blocklist = unlist(blocklist)
  # outer design matrix 'ODM' and matrix of outer weights 'W'
  ODM = matrix(0, mvs, lvs)
  for (j in 1:lvs) ODM[which(blocklist==j),j] = rep(1, blocks[j])
  W = ODM %*% diag(1/(apply(X%*%ODM, 2, sd)*sdv), lvs, lvs)
  w.old = rowSums(W)    
  w.dif <- itermax <- 1
  repeat 
  {            
    Y = X %*% W  # external estimation of LVs 'Y'
    Y = scale(Y) * sdv
    # matrix of inner weights 'e' 
    E <- switch(scheme, 
                "centroid" = sign(cor(Y) * (IDM + t(IDM))),
                "factor" = cor(Y) * (IDM + t(IDM)),
                "path" = get_path_scheme(IDM, Y))
    Z = Y %*% E  # internal estimation of LVs 'Z'
    Z = Z %*% diag(1/(apply(Z,2,sd)*sdv), lvs, lvs)  # scaling Z
    # computing outer weights 'w'
    for (j in 1:lvs)
    {
      X.blok = X[,which(blocklist==j)] 
      if (modes[j]=="A")# reflective way
        ODM[which(blocklist==j),j] <- (1/nrow(X)) * Z[,j] %*% X.blok
      if (modes[j]=="B")# formative way
        ODM[which(blocklist==j),j] <- solve.qr(qr(X.blok), Z[,j])
    }
    W = ODM
    w.new = rowSums(W)                
    w.dif = sum((abs(w.old) - abs(w.new))^2)  # difference of out.weights 
    if (w.dif < tol || itermax == iter) break
    w.old = w.new
    itermax = itermax + 1
  } # end repeat       
  W = ODM %*% diag(1/(apply(X %*% ODM, 2, sd)*sdv), lvs, lvs)
  w.new = rowSums(W)                
  names(w.new) = colnames(X)
  dimnames(W) = list(colnames(X), rownames(IDM))       
  res.ws = list(w.new, W, itermax)
  if (itermax == iter) res.ws = NULL
  return(res.ws)
}
