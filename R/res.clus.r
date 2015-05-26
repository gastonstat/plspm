#' @title Clustering on communality and structural residuals
#' 
#' @description
#' Computes communality and structural residuals from a global PLS-PM model 
#' and performs a Hierarchical Cluster Analysis on these residuals 
#' according to the REBUS algorithm.
#' 
#' @details
#' \code{res.clus()} comprises the second and third steps of 
#' the REBUS-PLS Algorithm. It computes communality and structural 
#' residuals. Then it performs a Hierarchical Cluster Analysis on these 
#' residuals (step three of REBUS-PLS Algorithm). As a result, this function 
#' directly provides a dendrogram obtained from a Hierarchical Cluster Analysis. 
#' 
#' @param pls Object of class \code{"plspm"}
#' @param Y Optional dataset (matrix or data frame) used when argument 
#' \code{dataset=NULL} inside \code{pls}.
#' @return An Object of class \code{"hclust"} containing the results of the
#' Hierarchical Cluster Analysis on the communality and structural residuals.
#' @references Esposito Vinzi V., Trinchera L., Squillacciotti S., 
#' and Tenenhaus M. (2008) REBUS-PLS: A Response-Based Procedure for 
#' detecting Unit Segments in PLS Path Modeling. \emph{Applied Stochastic Models 
#' in Business and Industry (ASMBI)}, \bold{24}, pp. 439-458. 
#' 
#' Trinchera, L. (2007) Unobserved Heterogeneity in Structural Equation Models: 
#' a new approach to latent class detection in PLS Path Modeling. 
#' \emph{Ph.D. Thesis}, University of Naples "Federico II", Naples, Italy.
#' @author Laura Trinchera, Gaston Sanchez
#' @seealso \code{\link{it.reb}}, \code{\link{plspm}}
#' @export 
#' @examples
#' \dontrun{
#'  ## example of rebus analysis with simulated data
#'    
#'  # load data
#'  data(simdata)
#'  
#'  # Calculate plspm
#'  sim_path = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
#'  dimnames(sim_path) = list(c("Price", "Quality", "Satisfaction"),
#'                             c("Price", "Quality", "Satisfaction"))
#'  sim_blocks = list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13)) 
#'  sim_modes = c("A", "A", "A")
#'  sim_global = plspm(simdata, sim_path, 
#'                     sim_blocks, modes=sim_modes)
#'  sim_global
#'    
#'  # Then compute cluster analysis on the residuals of global model
#'  sim_clus = res.clus(sim_global)
#'  }
#'
res.clus <- function(pls, Y = NULL)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm") 
    stop("\n'res.clus()' requires a 'plspm' object")
  # checking reflective modes
  if (any(pls$model$specs$modes != "A"))
    stop("\nSorry, REBUS only works for mode 'A'")
  # checking scaled data
  if (!pls$model$specs$scaled)
    stop("\nSorry, REBUS only works with scaled='TRUE'")
  # test availibility of dataset (either Y or pls$data)
  test_dataset(Y, pls$data, pls$model$gens$obs)
  
  # =======================================================
  # inputs setting
  # =======================================================
  IDM <- pls$model$IDM
  blocks <- pls$model$blocks
  blocklist = indexify(blocks)

  # data matrix DM
  if (!is.null(pls$data)) {
    DM = pls$data
    dataset = TRUE
  } else {         
    dataset = FALSE
    # building data matrix 'DM'
    DM = get_manifests(Y, blocks)
  }
  lvs = nrow(IDM)
  lvs.names = rownames(IDM)
  mvs = pls$model$gen$mvs
  # apply the selected scaling
  X = get_data_scaled(DM, TRUE)
    
  # =======================================================
  # computation of residuals
  # =======================================================
  Y.lvs <- pls$scores
  loads <- pls$outer_model$loading
  Path <- pls$path_coefs
  endo <- rowSums(IDM)
  endo[endo != 0] <- 1
  # matrices for storing outer and inner residuals
  outer_residuals = DM 
  inner_residuals = Y.lvs[,endo==1]
  # computation of outer residuals
  for (j in 1:lvs)
  {
    X.hat = Y.lvs[,j] %*% t(loads[blocklist==j])
    # outer residuals
    outer_residuals[,blocklist==j] = X[,blocklist==j] - X.hat
  }
  # computation of inner residuals
  # more than 1 endogenous LV
  if (sum(endo) != 1)
    Y.hat <- Y.lvs %*% t(Path[endo==1,])        
  # only 1 endogenous LV
  if (sum(endo) == 1)
    Y.hat = Y.lvs %*% Path[endo==1,]        
  # inner residuals
  inner_residuals = Y.lvs[,endo==1] - Y.hat
  
  # =======================================================
  # cluster analysis
  # =======================================================
  # hierarchical cluster analysis with Ward method using function "hcluster"
  res = cbind(outer_residuals, inner_residuals)    
  res.clus <- hcluster(res, method="euclidean", diag=FALSE, upper=FALSE,
                      link="ward", members=NULL, nbproc=2, doubleprecision=TRUE)
  # plot of the dendrogram
  plot(res.clus, main=c("REBUS", "Dendrogram of Outer and Inner Residuals"),
       hang=-1, cex.main=.9, cex.axis=.5, xlab="Hierarchical Clustering", 
       sub="Ward method", labels=FALSE)
  return(res.clus)
}
