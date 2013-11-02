#' @title Two Groups Comparison in PLS-PM
#' 
#' @description
#' Performs a group comparison test for comparing path coefficients 
#' between two groups. The null and alternative hypotheses to be tested are:
#' H0: path coefficients are not significantly different;
#' H1: path coefficients are significantly different
#' 
#' @details
#' \code{plspm.groups} performs a two groups comparison test in PLS-PM
#' for comparing path coefficients between two groups. Only two
#' methods are available: 1) bootstrap, and 2) permutation.
#' The bootstrap test is an adapted t-test based on bootstrap standard errors.
#' The permutation test is a randomization test which 
#' provides a non-parametric option.
#' 
#' When the object \code{pls} does not contain a data matrix 
#' (i.e. \code{pls$data=NULL}), the user must provide the data matrix 
#' or data frame in \code{Y}.
#' 
#' @param pls object of class \code{"plspm"}
#' @param group factor with 2 levels indicating the groups to be compared
#' @param Y optional dataset (matrix or data frame) used when argument
#' \code{dataset=NULL} inside \code{pls}.
#' @param method method to be used in the test. Possible values are 
#' \code{"bootstrap"} or \code{"permutation"}
#' @param reps integer indicating the number of either bootstrap resamples 
#' or number of permutations. If \code{NULL} then \code{reps}=100
#' @return An object of class \code{"plspm.groups"}
#' @return \item{test}{Table with the results of the applied test. 
#' Includes: path coefficients of the global model, path coeffs of group1, 
#' path coeffs of group2, (absolute) difference of path coeffs between groups, 
#' and the test results with the p-value.}
#' @return \item{global}{List with inner model results for the global model}
#' @return \item{group1}{List with inner model results for group1}
#' @return \item{group2}{List with inner model results for group2}
#' @author Gaston Sanchez
#' 
#' @references Chin, W.W. (2003) A permutation procedure for multi-group 
#' comparison of PLS models. In: Vilares M., Tenenhaus M., Coelho P., 
#' Esposito Vinzi V., Morineau A. (Eds.) \emph{PLS and Related Methods - 
#' Proceedings of the International Symposium PLS03.} Decisia, pp. 33-43.
#' 
#' Chin, W.W. (2000) Frequently Asked Questions, Partial Least Squares 
#' PLS-Graph. Available from: 
#' \url{http://disc-nt.cba.uh.edu/chin/plsfaq/multigroup.htm}
#' @seealso \code{\link{plspm}}
#' @export
#' @examples
#' 
#' \dontrun{
#'  ## example with customer satisfaction analysis
#'  ## group comparison based on the segmentation variable "gender"
#'  
#'  # load data satisfaction
#'  data(satisfaction)
#'  
#'  # define inner model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0) 
#'  LOY = c(1,0,0,0,1,0)
#'  sat_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'  
#'  # define outer model list
#'  sat_blocks = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'  
#'  # define vector of reflective modes
#'  sat_mod = rep("A", 6)
#'  
#'  # apply plspm
#'  satpls = plspm(satisfaction, sat_path, sat_blocks, 
#'                 modes = sat_mod, scaled = FALSE)
#'  
#'  # permutation test with 100 permutations
#'  group_perm = plspm.groups(satpls, satisfaction$gender, 
#'                            method="permutation", reps=100)
#'  group_perm
#'  }
#'
plspm.groups <-
function(pls, group, Y = NULL, method = "bootstrap", reps = NULL)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm") 
    stop("\n'pls' must be an object of class 'plspm'")
  g = group
  if (!is.factor(g)) 
    stop("\n'group' must be a factor")
  ng = nlevels(g)
  if (ng > 2) 
    stop("\n'group' must contain only 2 levels") 
  # test availibility of dataset (either Y or pls$data)
  test_dataset(Y, pls$data, pls$model$gens$obs)
  # check method
  if (!is.na(pmatch(method, "bootstrap"))) 
    method <- "bootstrap"
  METHODS <- c("bootstrap", "permutation")
  method <- pmatch(method, METHODS)
  if (is.na(method)) {
    warning("Invalid argument 'method'. Default 'method=bootstrap' is used.")   
    method <- "bootstrap"
  }
  # check number of replicates
  if (is.null(reps) | length(reps) > 1) reps = 100
  if (!is.numeric(reps) | floor(reps) <= 0) reps = 100
  
  # =======================================================
  # inputs setting
  # =======================================================  
  IDM = pls$model$IDM
  blocks = pls$model$blocks
  scaled = pls$model$specs$scaled
  blocklist  = indexify(blocks)
  
  # specifications
  specs = pls$model$specs
  
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
  X = get_data_scaled(DM, scaled)
  
  # =======================================================
  # Global model estimation
  # =======================================================  
  # object 'weights' contains outer w's, W, ODM, iter
  weights = get_weights(X, IDM, blocks, specs)
  ok_weights = test_null_weights(weights, specs)
  outer_weights = weights$w
  Y.lvs = get_scores(X, weights$W)
  
  # Path coefficients and total effects
  inner_results = get_paths(IDM, Y.lvs)
  innmod = inner_results[[1]]
  Path.global = inner_results[[2]]
  R2.global = inner_results[[3]]
  endo = rowSums(IDM)
  endo[endo != 0] = 1  # vector indicating endogenous LVs
  path.labs = NULL
  efs.labs = NULL
  for (j in 1:lvs)
    for (i in j:lvs)
      if (IDM[i,j] == 1) 
        path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
  path.global <- as.vector(Path.global[which(IDM==1)])
  names(path.global) <- path.labs
  
  # =======================================================
  # Group1 model estimation
  # =======================================================  
  g1.lab <- levels(g)[1]
  group1 <- which(g==levels(g)[1])
  ng1 <- length(group1)
  # apply the selected scaling
  X.g1 = get_data_scaled(DM[group1,], scaled)
  wgs.g1 = get_weights(X.g1, IDM, blocks, specs)
  ok_weights = test_null_weights(wgs.g1, specs)  
  Y1.lvs = get_scores(X.g1, wgs.g1$W)
  # Path coefficients 
  pathmod.g1 <- get_paths(IDM, Y1.lvs)
  innmod.g1 <- pathmod.g1[[1]]
  Path.g1 <- pathmod.g1[[2]]
  R2.g1 <- pathmod.g1[[3]]    
  path.g1 <- as.vector(Path.g1[which(IDM==1)])
  names(path.g1) <- path.labs
  
  # =======================================================
  # Group2 model estimation
  # =======================================================  
  g2.lab <- levels(g)[2]
  group2 <- which(g==levels(g)[2])
  ng2 <- length(group2)
  # apply the selected scaling
  X.g2 = get_data_scaled(DM[group2,], scaled)
  wgs.g2 = get_weights(X.g2, IDM, blocks, specs)
  ok_weights = test_null_weights(wgs.g2, specs)  
  Y2.lvs = get_scores(X.g2, wgs.g2$W)
  # Path coefficients 
  pathmod.g2 <- get_paths(IDM, Y2.lvs)
  innmod.g2 <- pathmod.g2[[1]]
  Path.g2 <- pathmod.g2[[2]]
  R2.g2 <- pathmod.g2[[3]]    
  path.g2 <- as.vector(Path.g2[which(IDM==1)])
  names(path.g2) <- path.labs
  
  # =======================================================
  # Group Comparison
  # =======================================================  
  dif.orig = abs(path.g1 - path.g2)
  nb = round(reps)
  
  if (method == 1)   # bootstrap
  {
    BG1 = matrix(0, nb, sum(IDM))
    BG2 = BG1
    for (i in 1:nb)
    {
      samg1 <- sample(group1, ng1, replace=TRUE) 
      samg2 <- sample(group2, ng2, replace=TRUE)
      # apply the selected scaling
      X.g1 = get_data_scaled(DM[samg1,], scaled)
      X.g2 = get_data_scaled(DM[samg2,], scaled)
      # outer weights
      wgs.g1 = get_weights(X.g1, IDM, blocks, specs)
      wgs.g2 = get_weights(X.g2, IDM, blocks, specs)
      if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
      if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples")
      # LVs scores
      Y1.lvs = get_scores(X.g1, wgs.g1$W)
      Y2.lvs = get_scores(X.g2, wgs.g2$W)
      # path coefficients
      pathmod.g1 = get_paths(IDM, Y1.lvs)
      paths.g1 = pathmod.g1[[2]]    
      pathmod.g2 = get_paths(IDM, Y2.lvs)
      paths.g2 = pathmod.g2[[2]]
      BG1[i,] = as.vector(paths.g1[which(IDM==1)])
      BG2[i,] = as.vector(paths.g2[which(IDM==1)])
    }    
    path.difs <- abs(apply(BG1,2,mean) - apply(BG2,2,mean))
    SE1 <- apply(BG1, 2, var)
    SE2 <- apply(BG2, 2, var)
    names(path.global) <- path.labs
    t.stat <- rep(NA, sum(IDM))
    k1 <- ((ng1-1)^2) / (ng1+ng2-2)
    k2 <- ((ng2-1)^2) / (ng1+ng2-2)
    k3 <- sqrt(1/ng1 + 1/ng2)
    for (i in 1:sum(IDM))         
      t.stat[i] <- path.difs[i] / (sqrt(k1*SE1[i]+k2*SE2[i]) * k3)        
    p.val <- pt(t.stat, ng1+ng2-2, lower.tail=FALSE)
    signi.path <- rep("no",length(p.val))
    signi.path[p.val<0.05] <- "yes"
    res.path <- round(cbind(path.global, path.g1, path.g2, dif.orig, 
                            t.stat, df=rep(ng1+ng2-2,sum(IDM)), p.val), 4)
    res <- data.frame(res.path, signi.path)
    colnames(res) <- c("global", paste(rep("group",2),levels(g),sep="."), 
                       "diff.abs", "t.stat", "deg.fr", "p.value", "sig.05") 
  } else
  {
    dif.perm <- matrix(0, nb, sum(IDM))
    for (i in 1:nb)# multigroup permutation
    {
      permu <- sample(1:(ng1+ng2), ng1+ng2)
      samg1 <- permu[1:ng1]
      samg2 <- permu[(ng1+1):(ng1+ng2)]
      # apply the selected scaling
      X.g1 = get_data_scaled(DM[samg1,], scaled)
      X.g2 = get_data_scaled(DM[samg2,], scaled)
      wgs.g1 <- get_weights(X.g1, IDM, blocks, specs)
      wgs.g2 <- get_weights(X.g2, IDM, blocks, specs)
      if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
      if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples") 
      # LVs scores
      Y1.lvs = get_scores(X.g1, wgs.g1$W)
      Y2.lvs = get_scores(X.g2, wgs.g2$W)
      # path coefficients
      pathmod.g1 = get_paths(IDM, Y1.lvs)
      paths.g1 = pathmod.g1[[2]]    
      pathmod.g2 = get_paths(IDM, Y2.lvs)
      paths.g2 = pathmod.g2[[2]]
      pp1 = as.vector(paths.g1[which(IDM==1)])
      pp2 = as.vector(paths.g2[which(IDM==1)])
      dif.perm[i,] = abs(pp1 - pp2)
    }   
    s.perm <- dif.orig 
    for (i in 1:sum(IDM))         
      s.perm[i] <- length(which(dif.orig[i] < dif.perm[,i])) + 1
    p.val <- (1/(nb+1)) * s.perm 
    signi.path = rep("no", length(p.val))
    signi.path[p.val < 0.05] = "yes"
    res.path = round(cbind(path.global, path.g1, path.g2, dif.orig, p.val), 4)
    res <- data.frame(res.path, signi.path)
    colnames(res) = c("global", paste(rep("group",2), levels(g), sep="."), 
                      "diff.abs", "p.value", "sig.05") 
  }
  
  # =======================================================
  # Results
  # =======================================================  
  met <- switch(method, "1"="bootstrap", "2"="permutation")
  settings <- c(scaled=scaled, scheme=pls$model$specs$scheme, method=met)
  res = list(test = res, 
             global = innmod, 
             group1 = innmod.g1, 
             group2 = innmod.g2, 
             settings = settings, 
             reps = reps)
  class(res) = "plspm.groups"
  return(res)
}


#' @S3method print plspm.groups
print.plspm.groups <- function(x,...)
{
  cat("GROUP COMPARISON IN PLS-PM FOR PATH COEFFICIENTS", "\n\n")
  cat("Scale of Data:      ", x$settings[[1]], "\n")
  cat("Weighting Scheme:   ", x$settings[[2]], "\n")
  cat("Selected method:    ", x$settings[[3]], "\n")
  cat("Num of replicates:  ", x$reps, "\n\n")
  cat("$test", "\n")
  print(x$test, print.gap=2)
  cat("\n")
  cat("Inner models in the following objects:", "\n")
  cat("$global ", "\n")
  cat("$group1 ", "\n")
  cat("$group2 ", "\n")
}
