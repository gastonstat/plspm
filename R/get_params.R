#' @title Check parameters for \code{plspm} and \code{plspm.fit}
#' 
#' @description
#' Internal function. \code{get_params} is called by \code{plspm}.
#'
#' @param x numeric matrix or data frame containing the manifest variables.
#' @param inner square (lower triangular) boolean matrix for inner model.
#' @param outer List of vectors with column indices from \code{x} indicating 
#' the sets of manifest variables asociated to the latent variables.
#' @param modes character vector indicating the type of measurement.
#' @param scheme string indicating the type of inner weighting scheme.
#' @param scaled logical indicating whether scaling data is performed.
#' @param boot.val logical indicating whether bootstrap validation is performed.
#' @param br integer indicating the number bootstrap resamples.
#' @param plsr logical indicating whether to use pls regression for path coefs.
#' @param tol decimal value indicating the tolerance criterion for covergence.
#' @param iter integer indicating the maximum number of iterations.
#' @param dataset logical indicating whether the data matrix should be retrieved.
#' @return list of validated parameters for \code{plspm} and \code{plspm.fit}.
#' @keywords internal
#' @export
get_params <-
function(x, inner, outer, modes=NULL, scheme="centroid", scaled=TRUE,
           boot.val=FALSE, br=NULL, plsr=FALSE, tol=0.00001, iter=100, dataset=TRUE)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (!is.matrix(x) && !is.data.frame(x))
    stop("Invalid object 'x'. Must be a numeric matrix or data frame.")
  if (is.null(rownames(x))) 
    rownames(x) <- 1:nrow(x)
  if (is.null(colnames(x))) 
    colnames(x) <- paste("MV", 1:ncol(x), sep="")
  if (!is.matrix(inner))
    stop("Invalid argument 'inner'. Must be a matrix.")
  if (nrow(inner)!=ncol(inner))
    stop("Invalid argument 'inner'. Must be a square matrix.")
  for (j in 1:ncol(inner))
    for (i in 1:nrow(inner)) {
      if (i<=j)
        if (inner[i,j]!=0) 
          stop("argument 'inner' must be a lower triangular matrix")
      if (length(intersect(inner[i,j], c(1,0)))==0)
        stop("elements in 'inner' must be '1' or '0'")
    }
  if (is.null(dimnames(inner)))
    lvs.names <- paste("LV", 1:ncol(inner), sep="")
  if (!is.null(rownames(inner)))
    lvs.names <- rownames(inner)
  if (!is.null(colnames(inner)))
    lvs.names <- colnames(inner)
  if (!is.list(outer))
    stop("Invalid argument 'outer'. Must be a list.")
  if (length(outer) != nrow(inner))
    stop("Number of rows of 'inner' does not coincide with length of 'outer'.")
  if (is.null(modes)) {
    modes <- rep("A",length(outer))
    warning("Argument 'modes' missing. Default reflective 'modes' is used.")
  }
  if (length(outer) != length(modes)) {
    warning("Warning: Invalid length of 'modes'. Default reflective 'modes' is used.")
    modes <- rep("A", length(outer))
  }
  for (i in 1:length(modes))
    if (modes[i]!="A" && modes[i]!="B") modes[i]<-"A"
  if (!is.na(pmatch(scheme, "centroid"))) 
    scheme <- "centroid"
  SCHEMES <- c("centroid", "factor", "path")
  scheme <- pmatch(scheme, SCHEMES)
  if (is.na(scheme)) {
    warning("Warning: Invalid argument 'scheme'. Default 'scheme=centroid' is used.")   
    scheme <- "centroid"
  }
  if (!is.logical(scaled)) {
    warning("Warning: Invalid argument 'scaled'. Default 'scaled=TRUE' is used.")
    scaled <- TRUE
  }
  if (!is.logical(boot.val)) {
    warning("Warning: Invalid argument 'boot.val'. No bootstrap validation is done.")
    boot.val <- FALSE
  }   
  if (boot.val) {
    if (!is.null(br)) {        
      if (mode(br)!="numeric" || length(br)!=1 || (br%%1)!=0 ||
        br<100 || br>1000) {
        warning("Warning: Invalid argument 'br'. Default 'br=100' is used.")   
        br <- 100
      } 
    } else
      br <- 100
  }
  if (!is.logical(plsr)) plsr<-FALSE
  if (mode(tol)!="numeric" || length(tol)!=1 || tol<=0 || tol>0.001) {
    warning("Warning: Invalid argument 'tol'. Default 'tol=0.00001' is used.")   
    tol <- 0.00001
  } 
  if (mode(iter)!="numeric" || length(iter)!=1 || iter<100) {
    warning("Warning: Invalid argument 'iter'. Default 'iter=100' is used.")   
    iter <- 100
  } 
  if (!is.logical(dataset)) 
    dataset <- TRUE
  
  # =======================================================
  # list with verified arguments
  # =======================================================
  res = list(
    x = x,
    inner = inner,
    outer = outer,
    modes = modes,
    scheme = scheme,
    SCHEMES = SCHEMES,
    scaled = scaled,
    boot.val = boot.val,
    br = br,
    plsr = plsr,
    tol = tol,
    iter = iter,
    dataset = dataset)
  res
}
