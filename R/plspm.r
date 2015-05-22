#' @title PLS-PM: Partial Least Squares Path Modeling
#'
#' @description
#' Estimate path models with latent variables by partial least squares approach
#' (for both metric and non-metric data)
#'
#' @details
#' The function \code{plspm} estimates a path model by partial least squares
#' approach providing the full set of results. \cr
#'
#' The argument \code{path_matrix} is a matrix of zeros and ones that indicates
#' the structural relationships between latent variables. \code{path_matrix} 
#' must be a lower triangular matrix; it contains a 1 when column \code{j}
#' affects row \code{i}, 0 otherwise. \cr
#'
#' @param Data matrix or data frame containing the manifest variables.
#' @param path_matrix A square (lower triangular) boolean matrix representing 
#' the inner model (i.e. the path relationships between latent variables).
#' @param blocks list of vectors with column indices or column names
#' from \code{Data} indicating the sets of manifest variables forming 
#' each block (i.e. which manifest variables correspond to each block).
#' @param scaling optional argument for runing the non-metric approach; 
#' it is a list of string vectors indicating the type of 
#' measurement scale for each manifest variable specified in \code{blocks}.
#' \code{scaling} must be specified when working with non-metric variables.
#' Possible values: \code{"num"} (linear transformation, 
#' suitable for numerical variables), \code{"raw"} (no transformation), 
#' \code{"nom"} (non-monotonic transformation, suitable for nominal variables), 
#' and \code{"ord"} (monotonic transformation, suitable for ordinal variables).
#' @param modes character vector indicating the type of measurement for each
#' block. Possible values are: \code{"A", "B", "newA", "PLScore", "PLScow"}. 
#' The length of \code{modes} must be equal to the length of \code{blocks}.
#' @param scheme string indicating the type of inner weighting
#' scheme. Possible values are \code{"centroid"}, \code{"factorial"}, or
#' \code{"path"}.
#' @param scaled whether manifest variables should be standardized. 
#' Only used when \code{scaling = NULL}. When (\code{TRUE}, data is 
#' scaled to standardized values (mean=0 and variance=1). 
#' The variance is calculated dividing by \code{N} instead of \code{N-1}).
#' @param tol decimal value indicating the tolerance criterion for the
#' iterations (\code{tol=0.000001}). Can be specified between 0 and 0.001.
#' @param maxiter integer indicating the maximum number of iterations
#' (\code{maxiter=100} by default). The minimum value of \code{maxiter} is 100.
#' @param plscomp optional vector indicating the number of PLS components
#' (for each block) to be used when handling non-metric data 
#' (only used if \code{scaling} is provided)
#' @param boot.val whether bootstrap validation should be performed. 
#' (\code{FALSE} by default). 
#' @param br number bootstrap resamples. Used only
#' when \code{boot.val=TRUE}. When \code{boot.val=TRUE}, the default number of 
#' re-samples is 100.
#' @param dataset whether the data matrix used in the computations should be
#' retrieved (\code{TRUE} by default).
#' @return An object of class \code{"plspm"}. 
#' @return \item{outer_model}{Results of the outer model. Includes:
#' outer weights, standardized loadings, communalities, and redundancies}
#' @return \item{inner_model}{Results of the inner (structural) model. 
#' Includes: path coeffs and R-squared for each endogenous latent variable}
#' @return \item{scores}{Matrix of latent variables used to estimate the inner
#' model. If \code{scaled=FALSE} then \code{scores} are latent variables
#' calculated with the original data (non-stardardized).}
#' @return \item{path_coefs}{Matrix of path coefficients 
#' (this matrix has a similar form as \code{path_matrix})}
#' @return \item{crossloadings}{Correlations between the latent variables 
#' and the manifest variables (also called crossloadings)}
#' @return \item{inner_summary}{Summarized results of the inner model. 
#' Includes: type of LV, type of measurement, number of indicators, R-squared,
#' average communality, average redundancy, and average variance
#' extracted}
#' @return \item{effects}{Path effects of the structural relationships. 
#' Includes: direct, indirect, and total effects}
#' @return \item{unidim}{Results for checking the unidimensionality of blocks
#' (These results are only meaningful for reflective blocks)}
#' @return \item{gof}{Goodness-of-Fit index}
#' @return \item{data}{Data matrix containing the manifest variables used in the
#' model. Only available when \code{dataset=TRUE}}
#' @return \item{boot}{List of bootstrapping results; only available 
#' when argument \code{boot.val=TRUE}}
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
#' 
#' Russolillo, G. (2012) Non-Metric Partial Least Squares. \emph{Electronic 
#' Journal of Statistics}, \bold{6}, pp. 1641-1669.
#' \url{http://projecteuclid.org/euclid.ejs/1348665231}
#' @seealso \code{\link{innerplot}}, \code{\link{outerplot}}, 
#' @export
#' @examples
#' \dontrun{
#' ## typical example of PLS-PM in customer satisfaction analysis
#' ## model with six LVs and reflective indicators
#' 
#' # load dataset satisfaction
#' data(satisfaction)
#' 
#' # path matrix
#' IMAG = c(0,0,0,0,0,0)
#' EXPE = c(1,0,0,0,0,0)
#' QUAL = c(0,1,0,0,0,0)
#' VAL = c(0,1,1,0,0,0)
#' SAT = c(1,1,1,1,0,0) 
#' LOY = c(1,0,0,0,1,0)
#' sat_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#' 
#' # plot diagram of path matrix
#' innerplot(sat_path)
#' 
#' # blocks of outer model
#' sat_blocks = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#' 
#' # vector of modes (reflective indicators)
#' sat_mod = rep("A", 6)
#' 
#' # apply plspm
#' satpls = plspm(satisfaction, sat_path, sat_blocks, modes = sat_mod, 
#'    scaled = FALSE)
#'    
#' # plot diagram of the inner model
#' innerplot(satpls)
#' 
#' # plot loadings
#' outerplot(satpls, what = "loadings")
#' 
#' # plot outer weights
#' outerplot(satpls, what = "weights")
#' }
plspm <-
function(Data, path_matrix, blocks, modes = NULL, scaling = NULL,  
         scheme = "centroid", scaled = TRUE, tol = 0.000001, maxiter = 100, 
         plscomp = NULL, boot.val = FALSE, br = NULL, 
         dataset = TRUE)
{
  # =======================================================================
  # checking arguments
  # =======================================================================
  valid = check_args(Data=Data, path_matrix=path_matrix, blocks=blocks, 
                     scaling=scaling, modes=modes, scheme=scheme, 
                     scaled=scaled, tol=tol, maxiter=maxiter, 
                     plscomp=plscomp, boot.val=boot.val, br=br, 
                     dataset=dataset)
  
  Data = valid$Data
  path_matrix = valid$path_matrix
  blocks = valid$blocks
  specs = valid$specs
  boot.val = valid$boot.val
  br = valid$br
  dataset = valid$dataset
  
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
  inner_results = get_paths(path_matrix, LV)
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
  crossloadings = data.frame(xloads, row.names=1:gens$mvs)
  crossloadings$name = factor(gens$mvs_names, levels = unique(gens$mvs_names))
  crossloadings$block = factor(rep(gens$lvs_names, block_sizes),
                               levels = gens$lvs_names)
  crossloadings = crossloadings[,c('name','block',colnames(xloads))]
  
  # outer model data frame
  outer_model = data.frame(
    name = factor(gens$mvs_names, levels = unique(gens$mvs_names)),
    block = factor(rep(gens$lvs_names, block_sizes),
                   levels = gens$lvs_names),
    weight = outer_weights, 
    loading = loadings, 
    communality = communality,
    redundancy = redundancy,
    row.names = 1:gens$mvs)
  
  # Unidimensionality
  unidim = get_unidim(DM = MV, blocks = blocks, modes = specs$modes)
  
  # Summary Inner model
  inner_summary = get_inner_summary(path_matrix, blocks, specs$modes,
                                    communality, redundancy, R2)
  
  # GoF Index
  gof = get_gof(communality, R2, blocks, path_matrix)
  
  # =======================================================================
  # Results
  # =======================================================================
  # deliver dataset?
  if (dataset) data = MV else data = NULL
  # deliver bootstrap validation results? 
  bootstrap = FALSE
  if (boot.val) 
  {
    if (nrow(X) <= 10) {
      warning("Bootstrapping stopped: very few cases.") 
    } else { 
      bootstrap = get_boots(MV, path_matrix, blocks, specs, br)
    }
  }
  
  # list with model specifications
  model = list(IDM=path_matrix, blocks=blocks, specs=specs,
               iter=weights$iter, boot.val=boot.val, br=br, gens=gens)
  
  ## output
  res = list(outer_model = outer_model, 
             inner_model = inner_model,
             path_coefs = Path, 
             scores = LV,
             crossloadings = crossloadings, 
             inner_summary = inner_summary, 
             effects = Path_effects,
             unidim = unidim, 
             gof = gof, 
             boot = bootstrap, 
             data = data,
             manifests = X,
             model = model)
  class(res) = "plspm"
  return(res)
}


#'@S3method print plspm
print.plspm <- function(x, ...)
{
  cat("Partial Least Squares Path Modeling (PLS-PM)", "\n")
  cat(rep("-", 45), sep="")
  cat("\n   NAME            ", "DESCRIPTION")  
  cat("\n1  $outer_model    ", "outer model")
  cat("\n2  $inner_model    ", "inner model")
  cat("\n3  $path_coefs     ", "path coefficients matrix")
  cat("\n4  $scores         ", "latent variable scores")
  if (!inherits(x, "plspm.fit"))
  {
    cat("\n5  $crossloadings  ", "cross-loadings")
    cat("\n6  $inner_summary  ", "summary inner model")
    cat("\n7  $effects        ", "total effects")
    cat("\n8  $unidim         ", "unidimensionality")
    cat("\n9  $gof            ", "goodness-of-fit")
    cat("\n10 $boot           ", "bootstrap results")
    cat("\n11 $data           ", "data matrix")
  }
  cat("\n")
  cat(rep("-", 45), sep="")
  cat("\nYou can also use the function 'summary'", "\n\n")    
  invisible(x)
}
