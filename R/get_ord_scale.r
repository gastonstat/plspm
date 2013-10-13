#' @title Non-Metric Ordinal Scale
#' 
#' @details
#' Internal function. \code{get_ord_scale} is called by \code{plspm}.
#'
#' @note
#' This function replaces the elements of x by the the means of y conditioned 
#' to the levels of x
#' 
#' @param y vector of values
#' @param x vector of the natural number series 1, 2, ..., p obtained 
#' by means of the function 'get_rank'
#' @param Xdummy dummy matrix corresponding to x (obtained with get_dummy)
#' @return scaled matrix
#' @keywords internal
#' @template internals
#' @export
get_ord_scale <- function(y, x, Xdummy) 
{
  n <- length(x)
  p <- max(x, na.rm = T)
  #	# ===========  build the (p x p) matrix of zeros ====	
  #	Xdummy<-matrix(0,n,p)
  #	# ===========  put the ones ============
  #	for (k in 1:p) {
  #		Xdummy[x == k,k] = 1
  #  	}
  #  	# ===========  if there are NA, add them ============
  #	if (any(is.na(x))) {
  #		Xdummy[which(rowSums(Xdummy) == 0),] <- NA
  #	}
  # =====  building an initial vector of scaling parameters ======
  quant <- (tapply(y, x, mean, na.rm=TRUE))
  
  # =====  searching for monotonically increasing quantifications ======
  quant_incr <- quant
  Xdummy_incr <- Xdummy
  repeat {
    ncol_Xdummy_Old <- ncol(Xdummy_incr)
    # ===== if the monotony is not respected, merge the columns =====	
    for (k in 1:(ncol(Xdummy_incr)-1)) {
      if (quant_incr[k] > quant_incr[k+1]) {
        Xdummy_incr[,k+1] <- Xdummy_incr[,k] + Xdummy_incr[,k+1]
        Xdummy_incr <- as.matrix(Xdummy_incr[,-k])
        quant_incr <- c()
        for (k in 1:ncol(Xdummy_incr)) {
          quant_incr[k] <- sum((Xdummy_incr[,k])*y,na.rm=T)
          quant_incr[k] <- quant_incr[k]/sum(Xdummy_incr[,k], na.rm=T)
        }
        break
      }
    }
    if (ncol(Xdummy_incr) == 1 || (ncol(Xdummy_incr) == ncol_Xdummy_Old)) {break}
  }
  x_quant_incr <- Xdummy_incr %*% quant_incr
  
  var_incr <- var(x_quant_incr, na.rm = TRUE)
  
  # =====  searching for monotonically decreasing quantifications ======
  quant_decr <- quant
  Xdummy_decr <- Xdummy
  repeat {
    ncol_Xdummy_Old<-ncol(Xdummy_decr)
    # ===== if the monotony is not respected, merge the columns =====		
    for (k in 1:(ncol(Xdummy_decr)-1)) {
      if (quant_decr[k] < quant_decr[k+1]) {
        Xdummy_decr[,k+1] <- Xdummy_decr[,k] + Xdummy_decr[,k+1]
        Xdummy_decr <- as.matrix(Xdummy_decr[,-k])
        quant_decr <- c()
        for (k in 1:ncol(Xdummy_decr)) {
          quant_decr[k] <- sum((Xdummy_decr[,k])*y,na.rm=T)
          quant_decr[k] <- quant_decr[k]/sum(Xdummy_decr[,k], na.rm=T)
        }
        break
      }
    }
    if (ncol(Xdummy_decr) == 1 || (ncol(Xdummy_decr) == ncol_Xdummy_Old)) {break}
  }
  x_quant_decr <- Xdummy_decr %*% quant_decr
  var_decr <- var(x_quant_decr, na.rm = TRUE)
  
  # =====  choosing between increasing and decreasing ======
  if (var_incr < var_decr) {
    Xdummy <- Xdummy_decr
    quant <- -(quant_decr)	
    x_quant <- -(x_quant_decr)	
  }
  else {
    Xdummy <- Xdummy_incr
    quant <- quant_incr	
    x_quant <- x_quant_incr	
  }
  
  x_quant
  # ===========  just in the case you need of them ============
  # ===========  eta2 = correlation rato  ============
  #eta2 <- var(x_quant)/var(y)
  #eta2<-var(x_quant, na.rm = T)/var(y, na.rm = T)	
  #list(Xdummy=Xdummy,x_quant=x_quant, eta2=eta2, quant=quant)
}
