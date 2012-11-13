#'@title PLS-PM for global and local models
#'
#'@description
#'Calculates PLS-PM for global and local models from a given partition
#'
#'@param pls An object of class \code{"plspm"}
#'@param y One object of the following classes: \code{"rebus"}, 
#'\code{"integer"}, or \code{"factor"}, that provides the class partitions.
#'@param Y Optional dataset (matrix or data frame) used when 
#'argument \code{dataset=NULL} inside \code{pls}.
#'
#'@details
#'\code{local.models} calculates PLS-PM for the global model 
#'(i.e. over all observations) as well as PLS-PM for local models 
#'(i.e. observations of different partitions).
#'
#'When \code{y} is an object of class \code{"rebus"}, \code{local.models} 
#'is applied to the classes obtained from the REBUS algorithm.
#'
#'When \code{y} is an \code{integer} vector or a \code{factor}, 
#'the values or levels are assumed to represent the group to which each 
#'observation belongs. In this case, the function \code{local.models} 
#'calculates PLS-PM for the global model, as well as PLS-PM for 
#'each group (local models).
#'
#'When the object \code{pls} does not contain a data matrix 
#'(i.e. \code{pls$data=NULL}), the user must provide the data matrix 
#'or data frame in \code{Y}. 
#'
#'The original parameters \code{modes}, \code{scheme}, \code{scaled}, 
#'\code{tol}, and \code{iter} from the object \code{pls} are taken.
#'
#'@return An object of class \code{"local.models"}, basically a list of 
#'length \code{k+1}, where \code{k} is the number of classes. 
#'@return \item{glob.model}{PLS-PM of the global model}
#'@return \item{loc.model.1}{PLS-PM of segment (class) 1}
#'@return \item{loc.model.2}{PLS-PM of segment (class) 2}
#'@return \item{loc.model.k}{PLS-PM of segment (class) k}
#'@note Each element of the list is an object of class \code{"plspm"}. 
#'Thus, in order to examine the results for each local model, 
#'it is necessary to use the \code{summary} function.
#'@author Laura Trinchera, Gaston Sanchez
#'@seealso \code{\link{rebus.pls}}
#'@export
#'@examples
#'  \dontrun{
#'    ## example of rebus analysis
#'    # load data
#'    data(sim.data)
#'    
#'    # compute GLOBAL model
#'    sim.inner = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
#'    dimnames(sim.inner) = list(c("Price", "Quality", "Satisfaction"),
#'                               c("Price","Quality", "Satisfaction"))
#'    sim.outer = list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13)) 
#'    sim.mod = c("A","A","A")  ## reflective indicators
#'    sim.global = plspm(sim.data, inner=sim.inner, 
#'                       outer=sim.outer, modes=sim.mod)
#'    sim.global
#'    
#'    # Then compute cluster on residual from global model
#'    sim.res.clus = res.clus(sim.global)
#'    
#'    # To conclude run iteration algorithm
#'    rebus.sim = it.reb(sim.global, sim.res.clus, nk=2, 
#'                        stop.crit = 0.005, iter.max = 100 )
#'    
#'    # Computation of local models 
#'    local.rebus = local.models(sim.global, rebus.sim)
#'    
#'    # Display plspm summary for first local model 
#'    summary(local.rebus$loc.model.1)
#'  }
#'
local.models <-
function(pls, y, Y=NULL)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls)!="plspm") 
    stop("\n'pls' must be an object of class 'plspm'")
  if (!is.element(class(y), c("rebus","integer","factor")))   
    stop("\n'y' must be of class 'rebus', 'integer' or 'factor'")
  if (class(y)=="rebus") {
    if (length(y$segments)!=nrow(pls$latents))
      stop("\n'pls' and 'y' are incompatible")
  } else {
    if (length(y)!=nrow(pls$latents))
      stop("\n'pls' and 'y' are incompatible")
  }
  if (!is.null(Y)) # if Y available
  {
    if (is.null(pls$data))
    {
      if (!is.matrix(Y) && !is.data.frame(Y))
        stop("\n'Y' must be a numeric matrix or data frame.")
      if (nrow(Y)!=nrow(pls$latents))
        stop("\n'pls' and 'Y' are incompatible. Different number of rows.")
    }
  } else { # if no Y
    if (is.null(pls$data)) 
      stop("\n'Y' is missing. No dataset available.")
  }
  
  # =======================================================
  # inputs settings
  # =======================================================
  IDM <- pls$model$IDM # Inner Design Matrix
  blocks <- pls$model$blocks # cardinality of blocks
  modes <- pls$model$modes # measurement modes    
  plsr <- FALSE 
  tol <- pls$model$tol
  iter <- pls$model$iter
  scheme <- pls$model$scheme
  scaled <- pls$model$scaled
  tol <- pls$model$tol
  iter <- pls$model$iter
  outer <- pls$model$outer
  blocklist <- outer
  for (k in 1:length(blocks))
    blocklist[[k]] <- rep(k,blocks[k])
  blocklist <- unlist(blocklist)
  # data matrix DM
  if (!is.null(pls$data)) {
    DM <- pls$data
    dataset <- TRUE
  } else {         
    dataset <- FALSE
    # building data matrix 'DM'
    DM <- matrix(NA, nrow(pls$latents), sum(blocks))
    for (k in 1:nrow(IDM))
      DM[,which(blocklist==k)] <- as.matrix(Y[,outer[[k]]])
    dimnames(DM) <- list(rownames(pls$latents), names(pls$out.weights))
  }
  lvs <- nrow(IDM)
  lvs.names <- rownames(IDM)
  mvs <- sum(blocks)
  endo <- rowSums(IDM)
  endo[endo!=0] <- 1  
  end.ind <- cumsum(blocks)
  ini.ind <- cumsum(blocks) - blocks + 1
  new.sets <- as.list(1:lvs)
  for (j in 1:lvs)
    new.sets[[j]] <- ini.ind[j]:end.ind[j]
  if (class(y)=="rebus") {
    segments <- as.factor(y$segments)
  } else {
    segments <- as.factor(y)
  }
  n.clus <- length(table(segments))
  
  # =======================================================
  # final models computation (global and local models)
  # =======================================================
  skem <- switch(scheme, "centroid"="centroid", "factor"="factor")
  final.mod <- as.list(1:(n.clus+1))# final plspm models
  for (k in 1:(n.clus+1))
  {
    if (k==1) {
      # global model
      X <- DM
      final.mod[[1]] = plspm(X, IDM, new.sets, modes, skem, scaled, 
                             tol=tol, iter=iter, dataset=dataset)
    } else
    {
      units.k <- which(segments==levels(segments)[k-1])
      # local models
      X.k <- DM[units.k,]
      final.mod[[k]] = plspm(X.k, IDM, new.sets, modes, skem, scaled, 
                             tol=tol, iter=iter, dataset=dataset)
    }
  }
  names(final.mod) = c("glob.model",paste(rep("loc.model",n.clus), 1:n.clus, sep="."))
  class(final.mod) = "local.models"
  return(final.mod)
}

