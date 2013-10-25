#' @title get_PLSRdoubleQ
#' 
#' @details
#' Internal function. \code{get_PLSRdoubleQ}
#'
#' @param Y numeric matrix
#' @param Yc numeric matrix
#' @param X numeric matrix
#' @param Xc numeric matrix
#' @param ncomp number of components
#' @keywords internal
#' @template internals
#' @export
get_PLSRdoubleQ <-
function (Y = NULL, Yc = NULL, X = NULL, Xc = NULL, ncomp) 
{
  if (!is.null(Y)) {
    n <- nrow(Y)
    rownamesY <- rownames(Y)
  }
  else {
    n <- nrow(Yc)
    rownamesY <- rownames(Yc)
  }
  ncolX <- 0
  if (!is.null(X)) 
    ncolX <- ncol(as.matrix(X))
  ncolXc <- 0
  if (!is.null(Xc)) 
    ncolXc <- ncol(as.matrix(Xc))
  if (!is.null(Xc) && !is.null(X)) 
    p <- ncolX + ncolXc
  if (!is.null(Xc) && is.null(X)) 
    p <- ncolXc
  if (!is.null(X) && is.null(Xc)) 
    p <- ncolX
  ncolY <- 0
  if (!is.null(Y)) 
    ncolY <- ncol(as.matrix(Y))
  ncolYc <- 0
  if (!is.null(Yc)) 
    ncolYc <- ncol(as.matrix(Yc))
  if (!is.null(Yc) && !is.null(Y)) 
    q <- ncolY + ncolYc
  if (!is.null(Yc) && is.null(Y)) 
    q <- ncolYc
  if (!is.null(Y) && is.null(Yc)) 
    q <- ncolY
  a <- ncomp
  Q <- matrix(, n, ncolXc)
  Qy <- matrix(, n, ncolYc)
  W <- matrix(, p, a)
  rownames(W) <- c(colnames(X), colnames(Xc))
  U <- matrix(, n, a)
  U <- matrix(c(rep(c(1, rep(0, (n - 1))), a)), n, a)
  rownames(U) <- rownamesY
  T <- matrix(, n, a)
  rownames(T) <- rownamesY
  C <- matrix(, q, a)
  rownames(C) <- c(colnames(Y), colnames(Yc))
  P <- matrix(, p, a)
  rownames(P) <- c(colnames(X), colnames(Xc))
  W_star <- matrix(, p, a)
  B <- matrix(, p, q)
  Pcorr <- matrix(, p, a)
  rownames(Pcorr) <- c(colnames(X), colnames(Xc))
  Ccorr <- matrix(, q, a)
  rownames(Ccorr) <- c(colnames(Y), colnames(Yc))
  Tcorr <- matrix(, n, a)
  rownames(Tcorr) <- rownamesY
  Xi <- X
  Yi <- Y
  Xarray <- array(, c(n, p, a))
  Yarray <- array(, c(n, q, a))
  for (i in 0:(a - 1)) {
    ncicli <- 0
    repeat {
      Ustart <- U[, i + 1]
      if (i == 0) {
        if (!is.null(Xc)) {
          Q <- dummy.G(U[, i + 1], Xc)$Quant
          Q <- scale(Q)
          colnames(Q) <- colnames(Xc)
          if (!is.null(X)) {
            Xi <- cbind(Xi[, 1:ncolX], Q)
          }
          else {
            Xi <- Q
          }
        }
      }
      W[, i + 1] <- as.matrix((t(Xi) %*% U[, i + 1])/as.numeric(t(U[, 
                                                                    i + 1]) %*% U[, i + 1]))
      W[, i + 1] <- W[, i + 1]/sqrt(as.numeric(t(W[, i + 
                                                     1]) %*% W[, i + 1]))
      T[, i + 1] <- (Xi %*% W[, i + 1])/as.numeric(t(W[, 
                                                       i + 1]) %*% W[, i + 1])
      if (i == 0) {
        if (!is.null(Yc)) {
          Qy <- dummy.G(T[, i + 1], Yc)$Quant
          Qy <- scale(Qy)
          if (!is.null(Y)) {
            Yi <- cbind(Yi[, 1:ncolY], Qy)
          }
          else {
            Yi <- Qy
          }
        }
      }
      C[, i + 1] <- (t(Yi) %*% T[, i + 1])/as.numeric(t(T[, 
                                                          i + 1]) %*% T[, i + 1])
      U[, i + 1] <- (Yi %*% C[, i + 1])/as.numeric(t(C[, 
                                                       i + 1]) %*% C[, i + 1])
      conv <- max(abs(Ustart - U[, i + 1]))
      ncicli <- ncicli + 1
      if (conv < 1e-07 | ncicli > 149) {
        break
      }
    }
    P[, i + 1] <- t(Xi) %*% T[, i + 1]/as.numeric(t(T[, i + 
                                                        1]) %*% T[, i + 1])
    Xi <- Xi - (T[, i + 1] %*% t(P[, i + 1]))
    Xarray[, , i + 1] <- Xi
    Yi <- Yi - (T[, i + 1] %*% t(C[, i + 1]))
    Yarray[, , i + 1] <- Yi
  }
  W_star <- W
  rownames(Q) <- rownamesY
  colnames(Q) <- colnames(Xc)
  rownames(B) <- c(colnames(X), colnames(Xc))
  colnames(B) <- c(colnames(Y), colnames(Yc))
  rownames(W_star) <- c(colnames(X), colnames(Xc))
  if (!is.null(Y)) {
    newY <- cbind(Y, Qy)
  }
  else {
    newY <- Qy
  }
  if (!is.null(X)) {
    newX <- cbind(X, Q)
  }
  else {
    newX <- Q
  }
  R2X <- 1 - (sum(apply(as.matrix(Xarray[, , a]), 2, var))/sum(apply(newX, 
                                                                     2, var)))
  R2Y <- 1 - (sum(apply(as.matrix(Yarray[, , a]), 2, var))/sum(apply(newY, 
                                                                     2, var)))
  if (a > 1) {
    W_star <- W %*% solve(t(P) %*% W)
    Pcorr <- P %*% (diag(apply(T, 2, sd)))
    Ccorr <- C %*% (diag(apply(T, 2, sd)))
    diag_matr <- diag(1/(apply(T, 2, sd) * sqrt(n - 1)))
    Tcorr <- T %*% (diag_matr)
    IDYarray <- array(, c(n, 2, q))
    lista_Ymedie <- list()
    for (j in 1:q) {
      IDYarray[, 1, j] <- (T %*% (C[j, ]))/sqrt(1 + (Ccorr[j, 
                                                           2]/Ccorr[j, 1])^2)
      IDYarray[, 2, j] <- IDYarray[, 1, j] * (Ccorr[j, 
                                                    2]/Ccorr[j, 1])
    }
    if (!is.null(Yc)) {
      for (k in 1:ncolYc) {
        matrice_Ymedie <- matrix(, max(as.matrix(Yc)[, 
                                                     k]), 2)
        matrice_Ymedie[, 1] <- as.vector(tapply(IDYarray[, 
                                                         1, k + ncolY], Yc[, k], mean, na.rm = T))
        matrice_Ymedie[, 2] <- as.vector(tapply(IDYarray[, 
                                                         2, k + ncolY], Yc[, k], mean, na.rm = T))
        lista_Ymedie[[k]] <- matrice_Ymedie
      }
    }
    IDarray <- array(, c(n, 2, p))
    lista_medie <- list()
    for (j in 1:p) {
      IDarray[, 1, j] <- (T %*% (P[j, ]))/sqrt(1 + (Pcorr[j, 
                                                          2]/Pcorr[j, 1])^2)
      IDarray[, 2, j] <- IDarray[, 1, j] * (Pcorr[j, 2]/Pcorr[j, 
                                                              1])
    }
    if (!is.null(Xc)) {
      for (k in 1:ncolXc) {
        matrice_medie <- matrix(, max(as.matrix(Xc)[, 
                                                    k]), 2)
        matrice_medie[, 1] <- as.vector(tapply(IDarray[, 
                                                       1, k + ncolX], Xc[, k], mean, na.rm = T))
        matrice_medie[, 2] <- as.vector(tapply(IDarray[, 
                                                       2, k + ncolX], Xc[, k], mean, na.rm = T))
        lista_medie[[k]] <- matrice_medie
      }
    }
    B <- W_star %*% t(C)
    VIP <- matrix(, p, 1)
    rownames(VIP) <- rownames(B)
    for (j in 1:p) {
      SumRdW2 <- 0
      SumRd <- 0
      for (h in 1:a) {
        Rd <- 0
        for (k in 1:q) {
          Rd <- Rd + (cor(newY[, k], T[, h])^2)
        }
        SumRdW2 <- SumRdW2 + ((Rd) * (as.numeric(W[j, 
                                                   h])^2))
        SumRd <- SumRd + (Rd)
      }
      VIP[j, 1] <- sqrt(p * SumRdW2/SumRd)
    }
    list(Q = Q, Qy = Qy, U = U, T = T, C = C, P = P, W = W, 
         B = B, W_star = W_star, Pcorr = Pcorr, Ccorr = Ccorr, 
         Tcorr = Tcorr, Xarray = Xarray, Yarray = Yarray, 
         IDarray = IDarray, lista_medie = lista_medie, IDYarray = IDYarray, 
         lista_Ymedie = lista_Ymedie, VIP = VIP, R2X = R2X, 
         R2Y = R2Y)
  }
  else {
    B <- W %*% t(C)
    list(Q = Q, Qy = Qy, U = U, T = T, C = C, P = P, W = W, 
         W_star = W_star, B = B, Xarray = Xarray, Yarray = Yarray, 
         R2X = R2X, R2Y = R2Y)
  }
}
