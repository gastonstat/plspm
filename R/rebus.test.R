#'@title Permutation Test for REBUS Multi-Group Comparison
#'
#'@description 
#'Performs permutation tests for comparing pairs of groups from a REBUS object.
#'
#'@details 
#'A permutation test on path coefficients, loadings, and GoF index 
#'is applied to the classes obtained from REBUS, by comparing two
#'classes at a time. That is to say, a permutation test is applied 
#'on pair of classes. The number of permutations in each test is 100. 
#'In turn, the number of classes handled by \code{rebus.test} is limited to 6. 
#'
#'When \code{pls$data=NULL} (there is no data matrix), the user
#'must provide the data matrix or data frame in \code{Y}.
#'
#'@param pls Object of class \code{"plspm"} returned by \code{\link{plspm}}
#'@param reb Object of class \code{"rebus"} returned by either 
#'\code{\link{rebus.pls}} or \code{\link{it.reb}}.
#'@param Y Optional dataset (matrix or data frame) used when argument 
#'\code{dataset=NULL} inside \code{pls}.
#'@return An object of class \code{"rebus.test"}, basically a list 
#'containing the results of each pair of compared classes. 
#'In turn, each element of the list is also a list with the results 
#'for the path coefficients, loadings, and GoF index.
#'@author Laura Trinchera, Gaston Sanchez
#'
#'@references Chin, W.W. (2003) A permutation procedure for multi-group 
#'comparison of PLS models. In: Vilares M., Tenenhaus M., Coelho P., 
#'Esposito Vinzi V., Morineau A. (Eds.) \emph{PLS and Related Methods - 
#'Proceedings of the International Symposium PLS03.} Decisia, pp. 33-43.
#'@seealso \code{\link{rebus.pls}}, \code{\link{local.models}}
#'@export
#'@examples
#'
#'  \dontrun{
#'  ## typical example of PLS-PM in customer satisfaction analysis
#'  ## model with six LVs and reflective indicators
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
#'  # Cluster analysis on residuals of global model
#'  sim_clus = res.clus(sim_global)
#'  
#'  # Iterative steps of REBUS algorithm on 2 classes
#'  rebus_sim = it.reb(sim_global, sim_clus, nk=2, 
#'                     stop.crit=0.005, iter.max=100)
#'                    
#'  # apply rebus.test
#'  sim_permu = rebus.test(sim_global, rebus_sim)
#'  
#'  # inspect sim.rebus
#'  sim_permu
#'  sim_permu$test_1_2
#'  
#'  # or equivalently
#'  sim_permu[[1]]
#'  }
#'
rebus.test <-
function(pls, reb, Y=NULL)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm") 
    stop("\n'pls' must be an object of class 'plspm'")
  # checking reflective modes
  if (any(pls$model$modes!="A"))
    stop("\nSorry, REBUS only works for reflective modes")
  # checking scaled data
  if (!pls$model$scaled)
    stop("\nSorry, REBUS only works with scaled='TRUE'")
  if (class(reb)!="rebus") 
    stop("\n'reb' must be an object of class 'rebus'")
  if (length(reb$segments)!=nrow(pls$data))
    stop("\n'pls' and 'reb' are incompatible")
  if (length(table(reb$segments)) > 6)
    stop("\nthe number of classes in 'rebus.test' is limited to 6")
  # if Y available
  if (!is.null(Y))
  {
    if (is.null(pls$data))
    {
      if (!is.matrix(Y) && !is.data.frame(Y))
        stop("\n'Y' must be a numeric matrix or data frame")
      if (nrow(Y) != nrow(pls$latents))
        stop("\n'pls' and 'Y' are incompatible. Different number of rows")
    }
  } else { # if no Y
    if (is.null(pls$data)) 
      stop("\n'Y' is missing. No dataset available.")
  }
  
  # =======================================================
  # Inputs settings
  # =======================================================
  IDM = pls$model$IDM # Inner Design Matrix
  blocks = pls$model$blocks # cardinality of blocks
  scheme = pls$model$scheme # inner weighting scheme
  modes = pls$model$modes # measurement modes
  scaled = pls$model$scaled # type of scaling
  plsr = pls$model$plsr # pls-regression
  tol = pls$model$tol # tolerance criterion
  iter = pls$model$iter # max num iterations
  outer = pls$model$outer
  blocklist = outer
  for (k in 1:length(blocks))
    blocklist[[k]] = rep(k, blocks[k])
  blocklist = unlist(blocklist)
  # data matrix DM
  if (!is.null(pls$data)) {
    DM = pls$data
    dataset = TRUE
  } else {         
    dataset = FALSE
    # building data matrix 'DM'
    DM = matrix(NA, nrow(pls$latents), sum(blocks))
    for (k in 1:nrow(IDM))
      DM[,which(blocklist==k)] = as.matrix(Y[,outer[[k]]])
    dimnames(DM) = list(rownames(pls$latents), names(pls$out.weights))
  }
  lvs = nrow(IDM)
  lvs.names = rownames(IDM)
  mvs = sum(blocks)
  # data scaling (standardized data)
  sd.X = sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
  X = scale(DM, scale=sd.X)
  n.clus = length(table(reb$segments))

  # =======================================================
  # multi-group comparison
  # =======================================================  
  ic = NULL
  ec = NULL
  for (i in 1:(n.clus-1))
  {
    ic = c(ic, rep(i, (n.clus-i)))
    ec = c(ec, seq((i+1), n.clus))
  }
  gp.index = cbind(ic, ec)
  gp.test = as.list(1:nrow(gp.index))
  for (i in 1:nrow(gp.index))
  { 
    a = which(reb$segments %in% gp.index[i,])
    g = as.factor(reb$segments[a])
    # apply locals test
    gp.test[[i]] = get_locals_test(DM[a,], pls, g)
  }
  # results
  names(gp.test) = paste(rep("test", nrow(gp.index)), 
                         gp.index[,1], gp.index[,2], 
                         sep = "_")
  class(gp.test) = "rebus.test"
  return(gp.test)
}

