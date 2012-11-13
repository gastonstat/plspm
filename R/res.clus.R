#'@title Cluster Analysis on communality and structural 
#'residuals of a PLS-PM object
#'
#'@description
#'Computes communality and structural residuals from the global model 
#'and performs a Hierarchical Cluster Analysis on these residuals 
#'according to the REBUS algorithm.
#'
#'@details
#'The function \code{res.clus} comprises the second and third steps of 
#'the REBUS-PLS Algorithm. It computes communality and structural 
#'residuals. Then it performs a Hierarchical Cluster Analysis on these 
#'residuals (step three of REBUS-PLS Algorithm). As a result, this function 
#'directly provides a dendrogram obtained from a Hierarchical Cluster Analysis. 
#'
#'@param pls Object of class \code{"plspm"}
#'@param Y Optional dataset (matrix or data frame) used when argument 
#'\code{dataset=NULL} inside \code{pls}.
#'@return An Object of class \code{"hclust"} containing the results of 
#'the Hierarchical Cluster Analysis on the communality and structural residuals.
#'@references Esposito Vinzi V., Trinchera L., Squillacciotti S., 
#'and Tenenhaus M. (2008) REBUS-PLS: A Response-Based Procedure for 
#'detecting Unit Segments in PLS Path Modeling. \emph{Applied Stochastic Models 
#'in Business and Industry (ASMBI)}, \bold{24}, pp. 439-458. 
#'
#'Trinchera, L. (2007) Unobserved Heterogeneity in Structural Equation Models: 
#'a new approach to latent class detection in PLS Path Modeling. 
#'\emph{Ph.D. Thesis}, University of Naples "Federico II", Naples, Italy.
#'@author Laura Trinchera, Gaston Sanchez
#'@seealso \code{\link{it.reb}}, \code{\link{plspm}}
#'@export 
#'@examples
#' \dontrun{
#'  ## example of rebus analysis with simulated data
#'    
#'  # load data
#'  data(simdata)
#'  
#'  # Calculate plspm
#'  sim_inner = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
#'  dimnames(sim_inner) = list(c("Price", "Quality", "Satisfaction"),
#'                             c("Price", "Quality", "Satisfaction"))
#'  sim_outer = list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13)) 
#'  sim_mod = c("A", "A", "A")  # reflective indicators
#'  sim_global = plspm(simdata, inner=sim_inner, 
#'                     outer=sim_outer, modes=sim_mod)
#'  sim_global
#'    
#'  # Then compute cluster analysis on the residuals of global model
#'  sim_clus = res.clus(sim_global)
#'  }
#'
res.clus <-
function(pls, Y=NULL)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls)!="plspm") 
    stop("\nAn object of class 'plspm' was expected")
  # checking reflective modes
  if (any(pls$model$modes != "A"))
    stop("\nSorry, REBUS only works for reflective modes")
  # checking scaled data
  if (!pls$model$scaled)
    stop("\nSorry, REBUS only works with scaled='TRUE'")
  # if Y available
  if (!is.null(Y))
  {
    if (is.null(pls$data))
    {
      if (!is.matrix(Y) && !is.data.frame(Y))
        stop("\n'Y'. Must be a numeric matrix or data frame.")
      if (nrow(Y)!=nrow(pls$latents))
        stop("\n'pls' and 'Y' are incompatible. Different number of rows.")
    }
  } else { 
    # if no Y
    if (is.null(pls$data)) 
      stop("\n'Y' is missing. No dataset available.")
  }
  
  # =======================================================
  # inputs setting
  # =======================================================
  IDM <- pls$model$IDM # Inner Design Matrix
  blocks <- pls$model$blocks # cardinality of blocks
  scheme <- pls$model$scheme # inner weighting scheme
  modes <- pls$model$modes # measurement modes
  scaled <- pls$model$scaled # type of scaling
  plsr <- pls$model$plsr # pls-regression
  tol <- pls$model$tol # tolerance criterion
  iter <- pls$model$iter # max num iterations
  outer <- pls$model$outer
  blocklist <- outer
  for (k in 1:length(blocks))
    blocklist[[k]] <- rep(k,blocks[k])
  blocklist <- unlist(blocklist)
  # data matrix DM
  if (!is.null(pls$data)) {
    DM <- pls$data
  } else {         
    # building data matrix 'DM'
    DM <- matrix(NA, nrow(pls$latents), sum(blocks))
    for (k in 1:nrow(IDM))
      DM[,which(blocklist==k)] <- as.matrix(Y[,outer[[k]]])
    dimnames(DM) <- list(rownames(pls$latents), names(pls$out.weights))
  }
  lvs <- nrow(IDM)
  lvs.names <- rownames(IDM)
  mvs <- sum(blocks)
  # data scaling (standardized data)
  sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
  X <- scale(DM, scale=sd.X)
  
  # =======================================================
  # computation of residuals
  # =======================================================
  Y.lvs <- pls$latents # recovering LV scores from pls
  loads <- pls$loadings # recovering loadings from pls
  PaCo <- pls$path.coefs # recovering path coeffs from pls
  endo <- rowSums(IDM)
  endo[endo!=0] <- 1  # indicator of endogenous LVs
  out.res <- DM # matrix for storing outer resids
  inn.res <- Y.lvs[,endo==1] # matrix for storing inner resids
  # computation of outer residuals
  for (j in 1:lvs)
  {
    q <- which(blocklist==j) 
    X.hat <- Y.lvs[,j] %*% t(loads[q])
    out.res[,q] <- X[,q] - X.hat # outer residuals
  }
  # computation of inner residuals
  # more than 1 endogenous LV
  if (sum(endo) != 1)
    Y.hat <- Y.lvs %*% t(PaCo[endo==1,])        
  # only 1 endogenous LV
  if (sum(endo) == 1)
    Y.hat <- Y.lvs %*% PaCo[endo==1,]        
  # inner residuals
  inn.res <- Y.lvs[,endo==1] - Y.hat
  
  # =======================================================
  # cluster analysis
  # =======================================================
  # hierarchical cluster analysis with Ward method using function "hcluster"
  res = cbind(out.res, inn.res)    
  res.clus = hcluster(res, method="euclidean", diag=FALSE, upper=FALSE,
                       link="ward", members=NULL, nbproc=2, doubleprecision=TRUE)
  # plot of the dendrogram
  plot(res.clus, main=c("REBUS", "Cluster Dendrogram of Outer and Inner Residuals"),
       hang=-1, cex.main=.9, cex.axis=.5, xlab="Hierarchical Clustering", 
       sub="Ward method", labels=FALSE)
  return(res.clus)
}
