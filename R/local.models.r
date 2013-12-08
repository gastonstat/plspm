#' @title PLS-PM for global and local models
#' 
#' @description
#' Calculates PLS-PM for global and local models from a given partition
#' 
#' @param pls An object of class \code{"plspm"}
#' @param y One object of the following classes: \code{"rebus"}, 
#' \code{"integer"}, or \code{"factor"}, that provides the class partitions.
#' @param Y Optional dataset (matrix or data frame) used when 
#' argument \code{dataset=NULL} inside \code{pls}.
#' 
#' @details
#' \code{local.models} calculates PLS-PM for the global model 
#' (i.e. over all observations) as well as PLS-PM for local models 
#' (i.e. observations of different partitions).
#' 
#' When \code{y} is an object of class \code{"rebus"}, \code{local.models} 
#' is applied to the classes obtained from the REBUS algorithm.
#' 
#' When \code{y} is an \code{integer} vector or a \code{factor}, 
#' the values or levels are assumed to represent the group to which each 
#' observation belongs. In this case, the function \code{local.models} 
#' calculates PLS-PM for the global model, as well as PLS-PM for 
#' each group (local models).
#' 
#' When the object \code{pls} does not contain a data matrix 
#' (i.e. \code{pls$data=NULL}), the user must provide the data matrix 
#' or data frame in \code{Y}. 
#' 
#' The original parameters \code{modes}, \code{scheme}, \code{scaled}, 
#' \code{tol}, and \code{iter} from the object \code{pls} are taken.
#' 
#' @return An object of class \code{"local.models"}, basically a list of 
#' length \code{k+1}, where \code{k} is the number of classes. 
#' @return \item{glob.model}{PLS-PM of the global model}
#' @return \item{loc.model.1}{PLS-PM of segment (class) 1}
#' @return \item{loc.model.2}{PLS-PM of segment (class) 2}
#' @return \item{loc.model.k}{PLS-PM of segment (class) k}
#' @note Each element of the list is an object of class \code{"plspm"}. 
#' Thus, in order to examine the results for each local model, 
#' it is necessary to use the \code{summary} function.
#' @author Laura Trinchera, Gaston Sanchez
#' @seealso \code{\link{rebus.pls}}
#' @export
#' @examples
#' \dontrun{
#' ## Example of REBUS PLS with simulated data
#' 
#' # load simdata
#' data("simdata", package='plspm')
#' 
#' # Calculate global plspm
#' sim_inner = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
#' dimnames(sim_inner) = list(c("Price", "Quality", "Satisfaction"),
#'                            c("Price", "Quality", "Satisfaction"))
#' sim_outer = list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13)) 
#' sim_mod = c("A", "A", "A")  # reflective indicators
#' sim_global = plspm(simdata, sim_inner, 
#'                    sim_outer, modes=sim_mod)
#' sim_global
#' 
#' ## Then compute cluster analysis on residuals of global model
#' sim_clus = res.clus(sim_global)
#' 
#' ## To complete REBUS, run iterative algorithm
#' rebus_sim = it.reb(sim_global, sim_clus, nk=2, 
#'                    stop.crit=0.005, iter.max=100)
#' 
#' ## You can also compute complete outputs 
#' ## for local models by running:
#' local_rebus = local.models(sim_global, rebus_sim)
#' 
#' # Display plspm summary for first local model
#' summary(local_rebus$loc.model.1)
#' }
local.models <- function(pls, y, Y = NULL)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm") 
    stop("\n'local.models()' requires an object of class 'plspm'")
  if (!is.element(class(y), c("rebus","integer","factor")))   
    stop("\n'y' must be of class 'rebus', 'integer' or 'factor'")
  if (class(y) == "rebus") {
    if (length(y$segments) != nrow(pls$scores))
      stop("\n'pls' and 'y' are incompatible")
  } else {
    if (length(y) != nrow(pls$scores))
      stop("\n'pls' and 'y' are incompatible")
  }
  # test availibility of dataset (either Y or pls$data)
  test_dataset(Y, pls$data, pls$model$gens$obs)
  
  # =======================================================
  # inputs settings
  # =======================================================
  IDM = pls$model$IDM
  blocks = pls$model$blocks
  scaling = pls$model$specs$scaling  
  modes = pls$model$specs$modes
  tol = pls$model$specs$tol
  maxiter = pls$model$specs$maxiter
  scheme = pls$model$specs$scheme
  scaled = pls$model$specs$scaled

  # data matrix DM
  if (!is.null(pls$data)) {
    DM = pls$data
    dataset = TRUE
  } else {         
    dataset = FALSE
    # building data matrix 'DM'
    DM = get_manifests(Y, blocks)
  }
  
  lvs = pls$model$gens$lvs
  lvs.names = pls$model$gens$lvs_names
  mvs = pls$model$gens$mvs
  endo = rowSums(IDM)
  endo[endo!=0] = 1
  aux_blocks = lengths(blocks)
  end.ind <- cumsum(aux_blocks)
  ini.ind <- cumsum(aux_blocks) - aux_blocks + 1
  new.sets <- as.list(1:lvs)
  for (j in 1:lvs)
    new.sets[[j]] <- ini.ind[j]:end.ind[j]
  if (class(y) == "rebus") {
    segments <- as.factor(y$segments)
  } else {
    segments <- as.factor(y)
  }
  n.clus <- length(table(segments))
  
  # =======================================================
  # final models computation (global and local models)
  # =======================================================
  # final plspm models
  final.mod <- as.list(1:(n.clus+1))
  for (k in 1:(n.clus+1))
  {
    if (k==1) {
      # global model
      final.mod[[1]] = plspm(DM, IDM, new.sets, scaling=scaling, modes=modes, 
                             scheme, scaled, tol, maxiter, dataset=dataset)
    } else
    {
      units.k <- which(segments == levels(segments)[k-1])
      # local models
      final.mod[[k]] = plspm(DM[units.k,], IDM, new.sets, scaling=scaling, 
            modes=modes, scheme, scaled, tol, maxiter, dataset=dataset)
    }
  }
  names(final.mod) = c("glob.model", 
                       paste(rep("loc.model", n.clus), 1:n.clus, sep="."))
  class(final.mod) = "local.models"
  return(final.mod)
}
