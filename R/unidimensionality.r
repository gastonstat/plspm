#' @title Unidimensionality of blocks
#' 
#' @description
#' Compute unidimensionality indices (a.k.a. Composite Reliability indices)
#'
#' @param Data matrix or data frame with variables
#' @param blocks optional list with vectors indicating the 
#' variables in each block
#' @return A data frame with the following columns:
#' @return \item{Block}{name of block}
#' @return \item{MVs}{number of manifest variables in each block}
#' @return \item{C.alpha}{Cronbach's alpha}
#' @return \item{DG.rho}{Dillon-Goldstein rho}
#' @return \item{eig.1st}{First eigenvalue}
#' @return \item{eig.2nd}{Second eigenvalue}
#' @author Gaston Sanchez
#' @seealso \code{\link{alpha}}, \code{\link{rho}}
#' @export
#' @examples
#'
#'  \dontrun{
#'  # load dataset satisfaction
#'  data(satisfaction)
#'
#'  # blocks Image and Expectations
#'  ima_expe = list(Image=1:5, Expec=6:10)
#'
#'  # compute unidimensionality indices
#'  unidim(satisfaction, ima_expe)
#'  }
#'
unidim <- function(Data, blocks = NULL)
{
  # check arguments
  DM = check_data(Data)
  if (!is.null(blocks)) {
    blocks = check_blocks(blocks, DM)
    DM = get_manifests(DM, blocks)
  } else {
    blocks = list(1:ncol(DM))
  }
  if (is.null(names(blocks)))
    names(blocks) = paste("block", 1:length(blocks), sep='')    
  
  # inputs setting
  lvs = length(blocks) 
  lvs_names = names(blocks)
  block_sizes = lengths(blocks)
  blocklist = indexify(blocks)
  correction = sqrt((nrow(DM)-1) / nrow(DM)) 
  # missing data flags
  missing_data = sapply(DM, is_missing)
  
  # initializing indices
  Alpha = rep(1, lvs)
  Rho = rep(1, lvs)
  eig.1st = rep(1, lvs)
  eig.2nd = rep(0, lvs)
  
  # calculate indices
  for (aux in 1:lvs)
  {
    if (any(missing_data[blocklist == aux]))
    {
      Alpha[aux] = NA
      Rho[aux] = NA
      eig.1st[aux] = NA
      eig.2nd[aux] = NA
    } else {
      if (block_sizes[aux] != 1) 
      { 
        # scaling block
        DM_block = DM[,blocklist == aux]
        stdev_X = apply(DM_block, 2, sd) * correction 
        X_scaled = scale(DM_block, scale=stdev_X)
        
        # PCA depending on block dimensions
        if (nrow(X_scaled) < ncol(X_scaled)) {   
          # more columns than rows
          X_pca = princomp(t(X_scaled)) 
          X_rho = t(X_scaled)
        } else {   
          # more rows than columns
          X_pca = princomp(X_scaled)
          X_rho = X_scaled
        }
        
        # indices
        Alpha[aux] = get_alpha(X_scaled)
        Rho[aux] = get_rho(X_rho, X_pca$scores[,1])
        eig.1st[aux] = X_pca$sdev[1]^2
        eig.2nd[aux] = X_pca$sdev[2]^2
      }
    }
  }
  
  # output
  data.frame(Block = lvs_names,
             MVs = block_sizes,
             C.alpha = Alpha, 
             DG.rho = Rho,
             eig.1st, 
             eig.2nd,
             row.names = NULL)
}


#' @title Cronbach's alpha
#' 
#' @description
#' Cronbach's alpha of a single block of variables
#'
#' @param X matrix representing one block of manifest variables 
#' @return Cronbach's alpha
#' @author Gaston Sanchez
#' @seealso \code{\link{rho}}, \code{\link{unidim}}
#' @export
#' @examples
#'
#'  \dontrun{
#'  # load dataset satisfaction
#'  data(satisfaction)
#'
#'  # block Image (first 5 columns of satisfaction)
#'  Image = satisfaction[,1:5]
#'
#'  # compute Cronbach's alpha for Image block
#'  alpha(Image)
#'  }
#'
alpha <- function(X)
{
  # checking input
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.numeric(X))
    stop("\n'alpha()' requires a numeric matrix")
  if (has_missing(X))
    stop("\n'X' contains missing values")
  
  # cronbach's alpha
  Alpha = 1
  if (ncol(X) > 1)
  {
    # correction factor
    correction = sqrt((nrow(X)-1) / nrow(X)) 
    # scaling data
    stdev_X = apply(X, 2, sd) * correction 
    X_scaled = scale(X, scale = stdev_X)
    # compute alpha
    Alpha = get_alpha(X_scaled)
  }
  
  # output
  Alpha
}


#' @title Dillon-Goldstein's rho
#' 
#' @description
#' Dillon-Goldstein's rho of a single block of variables
#'
#' @param X matrix representing one block of manifest variables 
#' @return Dillon-Goldstein's rho
#' @author Gaston Sanchez
#' @seealso \code{\link{alpha}}, \code{\link{unidim}}
#' @export
#' @examples
#'
#'  \dontrun{
#'  # load dataset satisfaction
#'  data(satisfaction)
#'
#'  # block Image (first 5 columns of satisfaction)
#'  Image = satisfaction[,1:5]
#'
#'  # compute Dillon-Goldstein's rho for Image block
#'  rho(Image)
#'  }
#'
rho <- function(X)
{
  # checking input
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.numeric(X))
    stop("\n'rho()' requires a numeric matrix")
  if (has_missing(X))
    stop("\n'X' contains missing values")
  
  # dillon-goldstein's rho
  Rho = 1
  if (ncol(X) > 1)
  {
    # correction factor
    correction = sqrt((nrow(X)-1) / nrow(X)) 
    # scaling data
    stdev_X = apply(X, 2, sd) * correction 
    X_scaled = scale(X, scale = stdev_X)
    
    # calculate rho
    if (nrow(X_scaled) < ncol(X_scaled)) {   
      # more columns than rows
      X_pca = princomp(t(X_scaled))
      Rho = get_rho(t(X_scaled), X_pca$scores[,1])
    } else {   
      # more rows than columns
      X_pca = princomp(X_scaled)
      Rho = get_rho(X_scaled, X_pca$scores[,1])
    }
  }
  
  # output
  Rho
}
