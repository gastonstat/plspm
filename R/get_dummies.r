#' @title Dummy matrices for categorical manifest variables
#' 
#' @details
#' Internal function. \code{get_dummies} is called by \code{plspm}.
#' 
#' @note If there are no categorical variables specified in specs,
#' \code{get_dummies} returns \code{NULL}
#' 
#' @param MV matrix or data frame from with manifest variables
#' @param specs list with algorithm specifications
#' @return dummies list with dummy matrices for categorical manifest variables
#' @keywords internal
#' @template internals
#' @export
get_dummies <- function(MV, specs)
{  
  # get metric
  metric = get_metric(specs$scaling)
  
  if (metric) {
    dummies = NULL
  } else {
    scaling = unlist(specs$scaling)
    nominal_ordinal = which(scaling %in% c("ord", "nom"))
    dummies = vector(mode="list", length(scaling))
    # only categorical mvs have an associated dummy matrix
    for (j in seq_along(nominal_ordinal)) {
      dummies[[nominal_ordinal[j]]] = get_dummy(MV[,nominal_ordinal[j]])
    }    
  }
  
  # output
  dummies
}


#' @title Non-Metric Dummy
#' 
#' @description
#' Internal function. \code{get_dummy} is called by \code{get_dummies}.
#' Transforms a vector of natural numbers from 1 to p into the
#' corresponding p x p dummy matrix
#' 
#' @param x a vector whose elements are non-negative integers from 0 to p
#' @return the dummy matrix
#' @keywords internal
#' @export
get_dummy <- function(x) 
{
  n = length(x)
  # p = max(x, na.rm = TRUE)
  # since 'x' could include zero, it's better to use the following:
  p = length(unique(x[!is.na(x)]))
  
  # build the (p x p) dummy matrix 
  Xdummy = matrix(0, n, p)
  for (k in 1:p) {
    Xdummy[x == k,k] = 1
  }
  # if there are NAs, add them
  if (any(is.na(x))) {
    Xdummy[which(rowSums(Xdummy) == 0),] <- NA
  }
  # output
  Xdummy
}


#' @title Dummy by Giorgio
#' 
#' @details
#' Internal function. \code{dummy.G} is called by \code{get_PLSRdoubleQ}.
#' 
#' @param Y a matrix
#' @param X a matrix
#' @return list with Xdummy, Quant, eta2
#' @keywords internal
#' @template internals
#' @export
dummy.G <- function (Y, X) 
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  p <- ncol(X)
  Ypred <- matrix(0, n, p)
  eta2 <- array(, p)
  for (k in 1:p) {
    Xtemp <- matrix(0, n, nlevels(as.factor(X[, k])))
    if (length(which(is.na(X[, k]))) > 0) 
      Xtemp[which(is.na(X[, k])), ] <- NA
    for (j in 1:nlevels(as.factor(X[, k]))) Xtemp[which(X[, 
                                                          k] == j), j] <- 1
    Ypred[, k] <- (Xtemp %*% (as.vector(tapply(Y, X[, k], 
                                               mean, na.rm = T))))
    eta2[k] <- var(Ypred[, k])/var(Y)
    if (k == 1) {
      Xdummy <- Xtemp
    }
    else {
      Xdummy <- cbind(Xdummy, Xtemp)
    }
  }
  list(Xdummy = Xdummy, Quant = Ypred, eta2 = eta2)
}
