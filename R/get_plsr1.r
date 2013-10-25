#' @title PLS regression for \code{plspm}
#' 
#' @details
#' Internal function. \code{get_plsr1} performs PLS-R1.
#' 
#' @param predictors matrix of predictors
#' @param response response variable
#' @param nc number of components
#' @param scaled logical indicating whether to scale the data
#' @return A list with pls regression results
#' @keywords internal
#' @template internals
#' @export
get_plsr1 <-
function(predictors, response, nc = NULL, scaled = TRUE)
{
  # ============ checking arguments ============
  X <- as.matrix(predictors)
  Y <- as.matrix(response)
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(nc))
    nc <- p
  # ============ setting inputs ==============
  if (scaled) Xx<-scale(X) else Xx<-scale(X,scale=F)
  if (scaled) Yy<-scale(Y) else Yy<-scale(Y,scale=F)
  X.old <- Xx
  Y.old <- Yy
  Th <- matrix(NA, n, nc)# matrix of X-scores
  Ph <- matrix(NA, p, nc)# matrix of X-loadings
  Wh <- matrix(NA, p, nc)# matrix of raw-weights
  Uh <- matrix(NA, n, nc)# matrix of Y-scores
  ch <- rep(NA, nc)# vector of y-loadings
  # ============ pls regression algorithm ==============
  for (h in 1:nc)
  {
    w.old <- t(X.old) %*% Y.old / sum(Y.old^2)
    w.new <- w.old / sqrt(sum(w.old^2)) # normalization
    t.new <- X.old %*% w.new
    p.new <- t(X.old) %*% t.new / sum(t.new^2) 
    c.new <- t(Y.old) %*% t.new / sum(t.new^2)
    u.new <- Y.old / as.vector(c.new)
    Y.old <- Y.old - t.new%*%c.new# deflate y.old
    X.old <- X.old - (t.new %*% t(p.new))# deflate X.old
    Th[,h] <- t.new
    Ph[,h] <- p.new
    Wh[,h] <- w.new
    Uh[,h] <- u.new
    ch[h] <- c.new
  }
  Ws <- Wh %*% solve(t(Ph)%*%Wh)# modified weights
  Bs <- as.vector(Ws %*% ch) # std beta coeffs    
  Br <- Bs * (rep(apply(Y, 2, sd),p)/apply(X,2,sd))   # beta coeffs
  cte <- as.vector(mean(response) - Br%*%apply(X,2,mean))# intercept
  y.hat <- X%*%Br+cte# y predicted
  resid <- as.vector(Y - y.hat)# residuals
  R2 <- as.vector(cor(Th, Yy))^2  # R2 coefficients    
  names(Br) <- colnames(X)
  names(resid) <- rownames(Y)
  names(y.hat) <- rownames(Y)
  names(R2) <- paste(rep("t",nc),1:nc,sep="")
  ## output
  list(coeffs = Br, 
       coef.std = Bs, 
       cte = cte, 
       R2 = R2[1:nc], 
       resid = resid, 
       y.pred = y.hat)
}
