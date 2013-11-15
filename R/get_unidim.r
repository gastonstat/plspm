#' @title Unidimensionality of reflective blocks
#' 
#' @details
#' Internal function. \code{get_unidim} is called by \code{plspm}.
#'
#' @param DM Data Matrix
#' @param blocks vector with numbers of variables per block
#' @param modes vector of modes
#' @return A data frame with the following columns:
#' @return \item{Type.measure}{Measurement mode}
#' @return \item{MVs}{number of manifest variables in each block}
#' @return \item{C.alpha}{Cronbach's alpha}
#' @return \item{DG.rho}{Dillon-Goldstein rho}
#' @return \item{eig.1st}{First eigenvalue}
#' @return \item{eig.2nd}{Second eigenvalue}
#' @keywords internal
#' @template internals
#' @export
get_unidim <- function(DM, blocks, modes)
{
  # inputs setting
  lvs = length(blocks) 
  lvs_names = names(blocks)
  blockinds = indexify(blocks)
  block_sizes = lengths(blocks)
  #  blocklist = unlist(lapply(block_sizes, function(x) rep(x, x)))
  obs = nrow(DM)
  sdvf = sqrt((nrow(DM)-1) / nrow(DM)) 
  # missing data flags
  missing_data = sapply(DM, is_missing)
  
  # Unidimensionality
  Alpha = rep(1, lvs)    # Cronbach's Alpha for each block
  Rho = rep(1, lvs)      # D.G. Rho for each block
  eig.1st = rep(1, lvs)  # first eigenvalue
  eig.2nd = rep(0, lvs)  # second eigenvalue

  # calculate indices
  for (aux in 1:lvs)
  {
    if (any(missing_data[blockinds == aux]))
    {
      Alpha[aux] = NA
      Rho[aux] = NA
      eig.1st[aux] = NA
      eig.2nd[aux] = NA
    } else {
      if (block_sizes[aux] != 1) 
      { 
        # scaling data
        DM.block = DM[,blockinds==aux]
        stdev.X = apply(DM.block, 2, sd) * sdvf 
        X_uni = scale(DM.block, scale=stdev.X)
        if (nrow(X_uni) < ncol(X_uni)) {   # more columns than rows
          acp = princomp(t(X_uni)) 
          X.rho = t(X_uni)
        } else {   # more rows than columns
          acp = princomp(X_uni)
          X.rho = X_uni
        }
        if (modes[aux] == "A") 
        {
          p = ncol(X_uni)
          # cronbach's alpha
          a.denom = var(rowSums(X_uni)) * sdvf^2
          a.numer = 2 * sum(cor(X_uni)[lower.tri(cor(X_uni))])
          alpha = (a.numer / a.denom) * (p / (p - 1))
          Alpha[aux] <- ifelse(alpha < 0, 0, alpha)
          # dillon-goldstein rho
          numer_rho <- colSums(cor(X.rho, acp$scores[,1]))^2
          denom_rho <- numer_rho + (p - colSums(cor(X.rho, acp$scores[,1])^2) )
          Rho[aux] <- numer_rho / denom_rho
        } else {  # modes[aux]=="B"
          Alpha[aux] = 0
          Rho[aux] = 0
        }
        eig.1st[aux] = acp$sdev[1]^2
        eig.2nd[aux] = acp$sdev[2]^2
      }      
    }
  }
  unidim = data.frame(Mode = modes, 
                      MVs = block_sizes,
                      C.alpha = Alpha, 
                      DG.rho = Rho,
                      eig.1st, 
                      eig.2nd)
  rownames(unidim) = lvs_names
  return(unidim)
}
