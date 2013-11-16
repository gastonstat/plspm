#' @title Basic results for Partial Least Squares Path Modeling
#' 
#' @description
#' Estimate path models with latent variables by partial least squares approach
#' without providing the full list of results as \code{plspm()}. This might be 
#' helpful when doing simulations, intensive computations, or when you don't 
#' want the whole enchilada.
#' 
#' @details
#' \code{plspm.fit} performs the basic PLS algorithm and provides limited
#' results (e.g. outer model, inner model, scores, and path coefficients). \cr
#' 
#' The argument \code{path_matrix} is a matrix of zeros and ones that indicates
#' the structural relationships between latent variables. \code{path_matrix} 
#' must be a lower triangular matrix; it contains a 1 when column \code{j}
#' affects row \code{i}, 0 otherwise. \cr
#' 
#' @param Data matrix or data frame containing the manifest variables.
#' @param path_matrix A square (lower triangular) boolean matrix representing 
#' the inner model (i.e. the path relationships betwenn latent variables).
#' @param blocks list of vectors with column indices or column names
#' from \code{Data} indicating the sets of manifest variables forming 
#' each block (i.e. which manifest variables correspond to each block).
#' @param scaling optional list of string vectors indicating the type of 
#' measurement scale for each manifest variable specified in \code{blocks}.
#' \code{scaling} must be specified when working with non-metric variables.
#' @param modes character vector indicating the type of measurement for each
#' block. Possible values are: \code{"A", "B", "newA", "PLScore", "PLScow"}. 
#' The length of \code{modes} must be equal to the length of \code{blocks}.
#' @param scheme string indicating the type of inner weighting
#' scheme. Possible values are \code{"centroid"}, \code{"factorial"}, or
#' \code{"path"}.
#' @param scaled whether manifest variables should be standardized. 
#' Only used when \code{scaling = NULL}. When (\code{TRUE} data is 
#' scaled to standardized values (mean=0 and variance=1). 
#' The variance is calculated dividing by \code{N} instead of \code{N-1}).
#' @param tol decimal value indicating the tolerance criterion for the
#' iterations (\code{tol=0.000001}). Can be specified between 0 and 0.001.
#' @param maxiter integer indicating the maximum number of iterations
#' (\code{maxiter=100} by default). The minimum value of \code{maxiter} is 100.
#' @param plscomp optional vector indicating the number of PLS components
#' (for each block) to be used when handling non-metric data 
#' (only used if \code{scaling} is provided)
#' @return An object of class \code{"plspm"}. 
#' @return \item{outer_model}{Results of the outer model. Includes:
#' outer weights, standardized loadings, communalities, and redundancies}
#' @return \item{inner_model}{Results of the inner (structural) model. 
#' Includes: path coeffs and R-squared for each endogenous latent variable}
#' @return \item{scores}{Matrix of latent variables used to estimate the inner
#' model. If \code{scaled=FALSE} then \code{scores} are latent variables
#' calculated with the original data (non-stardardized). If \code{scaled=TRUE}
#' then \code{scores} and \code{latents} have the same values}
#' @return \item{path_coefs}{Matrix of path coefficients 
#' (this matrix has a similar form as \code{path_matrix})}
#' @author Gaston Sanchez, Giorgio Russolillo
#'
#' @references Tenenhaus M., Esposito Vinzi V., Chatelin Y.M., and Lauro C.
#' (2005) PLS path modeling. \emph{Computational Statistics & Data Analysis},
#' \bold{48}, pp. 159-205.
#' 
#' Lohmoller J.-B. (1989) \emph{Latent variables path modeling with partial
#' least squares.} Heidelberg: Physica-Verlag.
#' 
#' Wold H. (1985) Partial Least Squares. In: Kotz, S., Johnson, N.L. (Eds.),
#' \emph{Encyclopedia of Statistical Sciences}, Vol. 6. Wiley, New York, pp.
#' 581-591.
#' 
#' Wold H. (1982) Soft modeling: the basic design and some extensions. In: K.G.
#' Joreskog & H. Wold (Eds.), \emph{Systems under indirect observations:
#' Causality, structure, prediction}, Part 2, pp. 1-54. Amsterdam: Holland.
#' @seealso \code{\link{innerplot}}, \code{\link{plot.plspm}}, 
#' @export
#' @examples
#'
#'  \dontrun{
#'  ## typical example of PLS-PM in customer satisfaction analysis
#'  ## model with six LVs and reflective indicators
#'
#'  # load dataset satisfaction
#'  data(satisfaction)
#'
#'  # inner model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0) 
#'  LOY = c(1,0,0,0,1,0)
#'  sat_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'
#'  # outer model list
#'  sat_blocks = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'
#'  # vector of reflective modes
#'  sat_modes = rep("A", 6)
#'
#'  # apply plspm.fit
#'  satpls = plspm.fit(satisfaction, sat_path, sat_blocks, sat_modes, 
#'      scaled=FALSE)
#'  
#'  # summary of results
#'  summary(satpls)
#'
#'  # default plot (inner model)
#'  plot(satpls)
#'  }
#'
plspm.fit <-
function(Data, path_matrix, blocks, modes = NULL, scaling = NULL, 
         scheme = "centroid", scaled = TRUE, tol = 0.000001,
         maxiter = 100, plscomp = NULL)
{
  # =======================================================================
  # checking arguments
  # =======================================================================
  valid = check_args(Data=Data, path_matrix=path_matrix, blocks=blocks, 
                     scaling=scaling, modes=modes, scheme=scheme, 
                     scaled=scaled, tol=tol, maxiter=maxiter, 
                     plscomp=plscomp, boot.val=FALSE, br=NULL, 
                     dataset=FALSE)
  
  Data = valid$Data
  path_matrix = valid$path_matrix
  blocks = valid$blocks
  specs = valid$specs
    
  # =======================================================================
  # Preparing data and blocks indexification
  # =======================================================================
  # building data matrix 'MV'
  MV = get_manifests(Data, blocks)
  check_MV = test_manifest_scaling(MV, specs$scaling)
  # generals about obs, mvs, lvs
  gens = get_generals(MV, path_matrix)
  # blocks indexing
  names(blocks) = gens$lvs_names
  block_sizes = lengths(blocks)
  blockinds = indexify(blocks)
  
  # transform to numeric if there are factors in MV
  if (test_factors(MV)) {
    numeric_levels = get_numerics(MV)
    MV = numeric_levels$MV
    categories = numeric_levels$categories
  }  
  # apply corresponding treatment (centering, reducing, ranking)
  X = get_treated_data(MV, specs)
  
  # =======================================================================
  # Outer weights and LV scores
  # =======================================================================
  metric = get_metric(specs$scaling)
  if (metric) {
    # object 'weights' contains outer w's, W, ODM, iter
    weights = get_weights(X, path_matrix, blocks, specs)
    ok_weights = test_null_weights(weights, specs)
    outer_weights = weights$w
    LV = get_scores(X, weights$W)
  } else {
    # object 'weights' contains outer w's, W, Y, QQ, ODM, iter
    weights = get_weights_nonmetric(X, path_matrix, blocks, specs)
    ok_weights = test_null_weights(weights, specs)
    outer_weights = weights$w
    LV = weights$Y
    X = weights$QQ  # quantified MVs
    colnames(X) = gens$mvs_names
  }
  
  # =======================================================================
  # Path coefficients and total effects
  # =======================================================================
  inner_results = get_paths(path_matrix, LV, FALSE)
  inner_model = inner_results[[1]]
  Path = inner_results[[2]]
  R2 = inner_results[[3]]
  Path_effects = get_effects(Path)
  
  # =======================================================================
  # Outer model: loadings, communalities, redundancy, crossloadings
  # =======================================================================
  xloads = cor(X, LV, use = 'pairwise.complete.obs')
  loadings = rowSums(xloads * weights$ODM)
  communality = loadings^2
  R2_aux = rowSums(weights$ODM %*% diag(R2, gens$lvs, gens$lvs))
  redundancy = communality * R2_aux
  
  # outer model data frame
  outer_model = data.frame(block = rep(gens$lvs_names, block_sizes),
                           weight = outer_weights, 
                           loading = loadings, 
                           communality = communality,
                           redundancy = redundancy,
                           row.names = gens$mvs_names)
  
  # =======================================================================
  # Results
  # =======================================================================
  # list with model specifications
  model = list(IDM=path_matrix, blocks=blocks, specs=specs, 
               iter=weights$iter, gens=gens)

  ## output
  res = list(outer_model = outer_model, 
             inner_model = inner_model,
             path_coefs = Path, 
             scores = LV,
             data = NULL, 
             model = model)
  class(res) = c("plspm.fit", "plspm")
  return(res)
}
