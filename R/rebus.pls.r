#' @title Response Based Unit Segmentation (REBUS)
#' 
#' @description
#' Performs all the steps of the REBUS-PLS algorithm. Starting from the 
#' global model, REBUS allows us to detect local models with better
#' performance.
#' 
#' @param pls Object of class \code{"plspm"}
#' @param Y Optional dataset (matrix or data frame) used when argument
#' \code{dataset=NULL} inside \code{pls}.
#' @param stop.crit Number indicating the stop criterion for the 
#' iterative algorithm. Use a threshold of less than 0.05\% of units 
#' changing class from one iteration to the other as stopping rule.
#' @param iter.max integer indicating the maximum number of iterations. 
#' @return An object of class \code{"rebus"}, basically a list with:
#' @return \item{loadings}{Matrix of standardized loadings 
#' (i.e. correlations with LVs.) for each local model.}
#' @return \item{path.coefs}{Matrix of path coefficients for each local model.}
#' @return \item{quality}{Matrix containing the average communalities, 
#' average redundancies, R2 values, and GoF values for each local model.}
#' @return \item{segments}{Vector defining for each unit the class membership.} 
#' @return \item{origdata.clas}{The numeric matrix with original data and 
#' with a new column defining class membership of each unit.}  
#' @author Laura Trinchera, Gaston Sanchez
#' @references Esposito Vinzi V., Trinchera L., Squillacciotti S., 
#' and Tenenhaus M. (2008) REBUS-PLS: A Response-Based Procedure for detecting 
#' Unit Segments in PLS Path Modeling. \emph{Applied Stochastic Models in 
#' Business and Industry (ASMBI)}, \bold{24}, pp. 439-458. 
#' 
#' Trinchera, L. (2007) Unobserved Heterogeneity in Structural Equation Models: 
#' a new approach to latent class detection in PLS Path Modeling. 
#' \emph{Ph.D. Thesis}, University of Naples "Federico II", Naples, Italy.
#' 
#' \url{http://www.fedoa.unina.it/2702/1/Trinchera_Statistica.pdf}
#' @seealso \code{\link{plspm}}, \code{\link{res.clus}}, 
#' \code{\link{it.reb}}, \code{\link{rebus.test}}, 
#' \code{\link{local.models}}
#' @export
#' @examples
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
#'  sim_global = plspm(simdata, sim_inner, 
#'                     sim_outer, modes=sim_mod)
#'  sim_global
#'  
#'  # run rebus.pls and choose the number of classes 
#'  # to be taken into account according to the displayed dendrogram.
#'  rebus_sim = rebus.pls(sim_global, stop.crit = 0.005, iter.max = 100)
#'  
#'  # You can also compute complete outputs for local models by running:
#'  local_rebus = local.models(sim_global, rebus_sim)
#'  }
#'
rebus.pls <-
function(pls, Y = NULL, stop.crit = 0.005, iter.max = 100)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm") 
    stop("\n'rebus.pls()' requires an object of class 'plspm'")
  # checking reflective modes
  if (any(pls$model$specs$modes != "A"))
    stop("\nREBUS only works with modes 'A'")
  # checking scaled data
  if (!pls$model$specs$scaled)
    stop("\nREBUS only works with scaled='TRUE'")
  # test availibility of dataset (either Y or pls$data)
  test_dataset(Y, pls$data, pls$model$gens$obs)
  if (mode(stop.crit) != "numeric" || length(stop.crit) != 1 || 
        stop.crit < 0 || stop.crit >= 1)
  {
    warning("Invalid stop criterion 'stop.crit'. Deafult value 0.005 is used")
    stop.crit = 0.005
  }
  if (mode(iter.max) != "numeric" || length(iter.max) != 1 || 
        iter.max <= 1 || (iter.max %% 1) != 0)
  {
    warning("Invalid 'iter.max'. Deafult value 100 is used")
    iter.max = 100
  }
  
  # =======================================================
  # perform REBUS algorithm
  # =======================================================
  resid = res.clus(pls, Y)
  print("Enter the number of classes (an integer > 1), and then press Enter:")
  nk <- scan(file="", n=1)
  if (mode(nk) != "numeric" || length(nk) != 1 || 
        nk <= 1 || (nk%%1) != 0)
    stop("\nInvalid number of classes. Must be an integer larger than 1")    
  resul = it.reb(pls, resid, nk, Y, stop.crit, iter.max)
  return(resul)
}
