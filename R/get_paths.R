#' @title Calculate path coefficients for \code{plspm}
#' 
#' @description
#' Internal function. \code{get_paths} is called by \code{plspm}.
#' 
#' @param IDM Inner Design Matrix
#' @param Y.lvs Matrix of latent variables
#' @param plsr logical to indicate path coefs by pls regression
#' @return list with inner mode, path coefs matrix, R2, and residuals
#' @keywords internal
#' @export
get_paths <-
  function(IDM, Y.lvs, plsr)
  {
    lvs.names <- colnames(IDM)
    endo = rowSums(IDM)
    endo[endo!=0] <- 1  # vector indicating endogenous LVs
    innmod <- as.list(1:sum(endo))
    Path <- IDM
    residuals <- as.list(1:sum(endo))
    R2 <- rep(0,nrow(IDM))
    for (aux in 1:sum(endo)) 
    {
      k1 <- which(endo==1)[aux]    # index for endo LV
      k2 <- which(IDM[k1,]==1)     # index for indep LVs
      if (length(k2)>1 & plsr) {               
        path.lm <- get_plsr1(Y.lvs[,k2], Y.lvs[,k1], nc=2)
        Path[k1,k2] <- path.lm$coeffs
        residuals[[aux]] <- path.lm$resid
        R2[k1] <- path.lm$R2[1]
        inn.val <- c(path.lm$R2[1], path.lm$cte, path.lm$coeffs)
        inn.lab <- c("R2", "Intercept", paste(rep("path_",length(k2)),names(k2),sep=""))
        names(inn.val) <- NULL
        innmod[[aux]] <- data.frame(concept=inn.lab, value=round(inn.val,4))
      }
      if (length(k2)==1 | !plsr) {
        path.lm <- summary(lm(Y.lvs[,k1] ~ Y.lvs[,k2]))
        Path[k1,k2] <- path.lm$coef[-1,1]
        residuals[[aux]] <- path.lm$residuals  
        R2[k1] <- path.lm$r.squared
        inn.val <- c(path.lm$r.squared, path.lm$coef[,1])
        inn.lab <- c("R2", "Intercept", paste(rep("path_",length(k2)),names(k2),sep=""))
        names(inn.val) <- NULL
        innmod[[aux]] <- data.frame(concept=inn.lab, value=round(inn.val,4))
      }
    }
    names(innmod) <- lvs.names[endo!=0]  
    names(R2) <- lvs.names
    res.paths <- list(innmod, Path, R2, residuals)
    return(res.paths)
  }
