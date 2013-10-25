#' @title Internal PLS regression (full data)
#' 
#' @details
#' Internal function. \code{get_PLSR}
#'
#' @param Y an (already centered) vector of order n
#' @param X an (already centered) n x p matrix
#' @param ncomp number of components
#' @keywords internal
#' @template internals
#' @export
# Output:
# \code{T} = X scores; \code{C} = Y weights; \code{P} = X loadings; 
# \code{A} = X weights; \code{B} = regression coefficients;
# \code{R2} = (2 x ncomp) matrix of the R2 of Y and X explained by models
# from 1 to \code{ncomp} components; \code{Pcorr} = cor(X,T);
# \code{Ccorr} = cor(Y,T); \code{Wstar} = X weights as \code{T=X%*%Wstar}
get_PLSR <- function(Y, X, ncomp) {
  Y = as.matrix(Y)
  X = as.matrix(X)
  colnamesX = colnames(X)
  n = nrow(Y)
  p = ncol(X)
  A = matrix(NA, p, ncomp)
  rownames(A) <- colnamesX
  colnames(A) <- c(1:ncomp)
  T <- matrix(NA, n, ncomp)
  rownames(T) <- rownames(X)
  colnames(T) <- c(1:ncomp)
  C <- matrix(NA, 1, ncomp)
  rownames(C) <- c("z")
  colnames(C) <- c(1:ncomp)
  P <- matrix(NA, p, ncomp)
  rownames(P) <- colnamesX
  colnames(P) <- c(1:ncomp)
  Wstar <- matrix(NA, p, ncomp)
  B <- matrix(,p,1)
  rownames(B) <- colnamesX
  colnames(B) <- c("z")
  R2 <- matrix(NA, 2, ncomp)
  varX <- sum(apply(X,2,var))
  varY <- var(Y)
  for (k in 1:ncomp) {
    if (k == 1) {
      A[,k] <- t(X) %*% Y
      A[,k] <- A[,k] / sqrt(sum(A[,k]^2))
      T[,k] <- X %*% A[,k]
      C[1,k] <- t(Y) %*% T[,k]/(sum(T[,k]^2))
      P[,k] <- t(X) %*% T[,k]/(sum(T[,k]^2))
      Xres <- X - T[,k] %*% t(P[,k])		
      if (ncomp == 1) {
        B <- A[,k] * C[1,k]
      }
      Wstar <- as.matrix(A[,k])
      R2[1,k] <- 1 - (sum(apply(Xres,2,var)) / varX)
      R2[2,k] <- 1 - (var(Y-(T[,k]%*%t(C[,k]))) / varY)
    }
    else {
      A[,k] <- t(Xres) %*% Y
      A[,k] <- A[,k] / sqrt(sum(A[,k]^2))
      T[,k] <- Xres %*% A[,k]
      C[1,k] <- t(Y) %*% T[,k]/(sum(T[,k]^2))
      P[,k] <- t(Xres) %*% T[,k]/(sum(T[,k]^2))
      Xres <- Xres - T[,k] %*% t(P[,k])
      R2[1,k] <- 1 - (sum(apply(Xres,2,var)) / varX)
      R2[2,k] <- 1 - (var(Y-(T[,1:k] %*% as.matrix(C[,1:k]))) / varY)
    }
  }
  Wstar <- A %*% solve(t(P) %*% A)	
  B <- Wstar %*% t(C)
  Pcorr <- as.matrix(cor(X,T))
  Ccorr <- as.matrix(cor(Y,T))
  rownames(R2) <- c("X","z")
  colnames(R2) <- c(1:ncomp)
  rownames(Pcorr) <- colnamesX
  colnames(Pcorr) <- c(1:ncomp) 
  rownames(Ccorr) <- c("z")
  colnames(Ccorr) <- c(1:ncomp)
  # output
  list(T = T, 
       C = C, 
       P = P, 
       A = A, 
       B = B, 
       R2cum = R2, 
       Pcorr = Pcorr, 
       Ccorr = Ccorr, 
       Wstar = Wstar)
}
