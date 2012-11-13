#' @title Unidimensionality of reflective blocks
#' 
#' @description
#' Internal this function. \code{get_unidim} is called by \code{plspm}
#'
#' @param x numeric matrix or data frame containing the manifest variables.
#' @param outer list of vectors with column indices from \code{x} indicating 
#' the sets of manifest variables asociated to the latent variables
#' @param modes character vector indicating the type of measurement
#' @param DM Data Matrix
#' @param blocks vector with numbers of variables per block
#' @param check logical to check \code{x}, \code{outer} and \code{modes}.
#' @return A data frame with the following columns:
#' @return \item{Type.measure}{Measurement mode}
#' @return \item{MVs}{number of manifest variables in each block}
#' @return \item{C.alpha}{Cronbach's alpha}
#' @return \item{DG.rho}{Dillon-Goldstein rho}
#' @return \item{eig.1st}{First eigenvalue}
#' @return \item{eig.2nd}{Second eigenvalue}
#' @keywords internal
#' @export
get_unidim <-
function(x=NULL, outer=NULL, modes=NULL, DM=NULL, blocks=NULL, check=TRUE)
{
  # =======================================================
  # checking first three arguments when called by users
  # =======================================================
  if (check)
  {
    # check x
    if (!is.matrix(x) && !is.data.frame(x))
      stop("Invalid object 'x'. Must be a numeric matrix or data frame.")
    if (is.null(rownames(x))) 
      rownames(x) = 1:nrow(x)
    if (is.null(colnames(x))) 
      colnames(x) = paste("MV", 1:ncol(x), sep="")
    # check outer
    if (!is.list(outer))
      stop("Invalid argument 'outer'. Must be a list.")
    if (is.null(modes)) {
      modes = rep("A", length(outer))
      warning("Argument 'modes' missing. Default reflective 'modes' is used.")
    }
    if (length(outer) != length(modes)) {
      warning("Warning: Invalid length of 'modes'. Reflective 'modes' is used.")
      modes = rep("A", length(outer))
    }
    for (i in 1:length(modes))
      if (modes[i] != "A" && modes[i] != "B") modes[i] = "A"
    # get blocks
    blocks = unlist(lapply(outer, length))
    mvs = sum(blocks)
    names(blocks) = paste("Block", 1:length(blocks), sep="")
  }
  
  # =======================================================
  # inputs setting
  # =======================================================  
  lvs = length(blocks) 
  lvs.names = names(blocks)
  blocklist = as.list(1:lvs)
  for (j in 1:lvs)
    blocklist[[j]] = rep(j,blocks[j])
  blocklist = unlist(blocklist)
  Mode <- modes
  Mode[modes=="A"] = "Reflective"
  Mode[modes=="B"] = "Formative"   
  obs = nrow(DM)
  sdvf = sqrt((nrow(DM)-1) / nrow(DM)) 
  # Unidimensionality
  Alpha = rep(1, lvs)   # Cronbach's Alpha for each block
  Rho = rep(1, lvs)     # D.G. Rho for each block
  eig.1st = rep(1, lvs) # first eigenvalue
  eig.2nd = rep(0, lvs) # second eigenvalue
  for (aux in 1:lvs) 
  {      
    if (blocks[aux] != 1) 
    { 
      # scaling data
      DM.block = DM[,which(blocklist==aux)]
      stdev.X = apply(DM.block, 2, sd) * sdvf 
      X.uni = scale(DM.block, scale=stdev.X)
      if (nrow(X.uni)<ncol(X.uni)) {   # more columns than rows
        acp = princomp(t(X.uni)) 
        X.rho = t(X.uni)
      } else {   # more rows than columns
        acp = princomp(X.uni)
        X.rho = X.uni
      }
      if (modes[aux]=="A") 
      {
        p = ncol(X.uni)
        # cronbach's alpha
        a.denom <- var(rowSums(X.uni)) * sdvf^2
        a.numer <- 2*sum(cor(X.uni)[lower.tri(cor(X.uni))])
        alpha <- (a.numer / a.denom) * (p/(p-1))
        Alpha[aux] <- ifelse(alpha<0,0,alpha)
        # dillon-goldstein rho
        numer.rho <- colSums(cor(X.rho, acp$scores[,1]))^2
        denom.rho <- numer.rho + (p - colSums(cor(X.rho, acp$scores[,1])^2) )
        Rho[aux] <- numer.rho / denom.rho
      } else {  # modes[aux]=="B"
        Alpha[aux] = 0
        Rho[aux] = 0
      }
      eig.1st[aux] = acp$sdev[1]^2
      eig.2nd[aux] = acp$sdev[2]^2
    }
  }
  unidim = data.frame(Type.measure = Mode, 
                      MVs = blocks,
                      C.alpha = Alpha, 
                      DG.rho = Rho,
                      eig.1st, 
                      eig.2nd)
  rownames(unidim) = lvs.names
  return(unidim)
}