#' @title Internal PLS regression with missing values
#' 
#' @details
#' Internal function. \code{get_PLSR_NA}
#'
#' @param Y an (already centered) vector of order n; \code{NA} not allowed
#' @param X an (already centered) n x p matrix; \code{NA} allowed but
#' with at least one observation for every row and column
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
get_PLSR_NA <- function(Y, X, ncomp) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  X_avail <- 1-is.na(X)
  colnamesX <- colnames(X)
  n <- nrow(Y)
  p <- ncol(X)
  A <- matrix(,p,ncomp)
  rownames(A) <- colnamesX
  colnames(A) <- c(1:ncomp)
  T <- matrix(,n,ncomp)
  rownames(T) <- rownames(X)
  colnames(T) <- c(1:ncomp)
  C <- matrix(,1,ncomp)
  rownames(C) <- c("z")
  colnames(C) <- c(1:ncomp)
  P <- matrix(,p,ncomp)
  rownames(P) <- colnamesX
  colnames(P) <- c(1:ncomp)
  Wstar <- matrix(,p,ncomp)
  B <- matrix(,p,1)
  rownames(B) <- colnamesX
  colnames(B) <- c("z")
  R2 <- matrix(,2,ncomp)
  varX <- sum(apply(X,2,function(x){var(x,na.rm=TRUE)}))
  varY <- var(Y)
  for (k in 1:ncomp) {
    if (k == 1) {
      A[,k] <- colSums(X*Y[,1], na.rm = TRUE)
      A[,k] <- A[,k]/colSums((X_avail*Y[,1])^2) 
      A[,k] <- A[,k]/sqrt(sum(A[,k]^2))
      T[,k] <- colSums(t(X)*A[,k], na.rm=TRUE) 
      T[,k] <- T[,k]/colSums((t(X_avail)*A[,k])^2)
      C[1,k] <- sum(Y[,1]*T[,k], na.rm=TRUE)/(sum(T[,k]^2, na.rm=TRUE))
      P[,k] <- colSums(X*T[,k],na.rm=TRUE)
      P[,k] <- P[,k]/colSums((X_avail*T[,k])^2) 
      Xres <- X - T[,k]%*%t(P[,k])
      Yres <- Y - T[,k]*C[1,k]		
      if (ncomp == 1) {
        B <- A[,k]*C[1,k]
      }
      Wstar <- as.matrix(A[,k])
      R2[1,k] <-1 - (sum(apply(Xres,2,function(x){var(x,na.rm=TRUE)})) / varX)
      R2[2,k] <- 1 - (var(Yres)/varY)
    }
    else {
      A[,k] <-colSums(Xres*Yres[,1],na.rm=TRUE)
      A[,k] <- A[,k]/colSums((X_avail*Yres[,1])^2) 
      A[,k] <- A[,k]/sqrt(sum(A[,k]^2))
      T[,k] <- colSums(t(Xres)*A[,k], na.rm=TRUE) 
      T[,k] <- T[,k]/colSums((t(X_avail)*A[,k])^2)
      C[1,k] <- sum(Yres[,1]*T[,k], na.rm=TRUE)/(sum(T[,k]^2, na.rm=TRUE))
      P[,k] <- colSums(Xres*T[,k],na.rm=TRUE)  
      P[,k] <- P[,k]/colSums((X_avail*T[,k])^2)
      Xres <- Xres - T[,k]%*%t(P[,k])
      Yres <- Yres - T[,k]*C[1,k]
      R2[1,k] <-1 - (sum(apply(Xres,2,function(x){var(x,na.rm=TRUE)})) / varX)
      R2[2,k] <- 1 - (var(Yres)/varY)
    }
  }
  uppertri_PA <- ((t(P)%*%A)*upper.tri(diag(ncomp)))+diag(ncomp)
  Wstar <- A%*%solve(uppertri_PA)	
  B <- Wstar%*%t(C)
  Pcorr <- as.matrix(cor(X,T,use="pairwise.complete.obs"))
  Ccorr <- as.matrix(cor(Y,T))
  rownames(Wstar) <- colnamesX
  colnames(Wstar) <- c(1:ncomp)
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
       Wstar = Wstar, 
       Xres = Xres, 
       Yres = Yres)
}
