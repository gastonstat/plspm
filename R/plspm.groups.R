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
#' @return \item{global}{List with results of the inner model for the global model}
#' @return \item{group1}{List with results of the inner model for group1}
#' @return \item{group2}{List with results of the inner model for group2}
#' @author Gaston Sanchez
#' 
#' @references Chin, W.W. (2003) A permutation procedure for multi-group comparison
#' of PLS models. In: Vilares M., Tenenhaus M., Coelho P., Esposito Vinzi V., 
#' Morineau A. (Eds.) \emph{PLS and Related Methods - Proceedings of the International 
#' Symposium PLS03.} Decisia, pp. 33-43.
#' 
#' Chin, W.W. (2000) Frequently Asked Questions, Partial Least Squares PLS-Graph. 
#' Available from: \url{http://disc-nt.cba.uh.edu/chin/plsfaq/multigroup.htm}
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
#'  sat.inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'  
#'  # define outer model list
#'  sat.outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'  
#'  # define vector of reflective modes
#'  sat.mod = rep("A", 6)
#'  
#'  # apply plspm
#'  my_pls = plspm(satisfaction, sat.inner, sat.outer, sat.mod, scheme="factor", 
#'               scaled=FALSE)
#'               
#'  # permutation test with 100 permutations
#'  group_comp = plspm.groups(my_pls, satisfaction$gender, 
#'                            method="permutation", reps=100)
#'  group_comp
#'  }
#'
plspm.groups <-
function(pls, group, Y = NULL, method = "bootstrap", reps = NULL)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm") 
    stop("'pls' must be an object of class 'plspm'")
  g = group
  if (!is.factor(g)) stop("'group' must be a factor")
  ng = nlevels(g)
  if (ng > 2) stop("'group' must contain only 2 levels") 
  if (!is.null(Y)) # if Y available
  {
    if (is.null(pls$data))
    {
      if (!is.matrix(Y) && !is.data.frame(Y))
        stop("Invalid object 'Y'. Must be a numeric matrix or data frame.")
      if (nrow(Y)!=nrow(pls$latents))
        stop("Argument 'pls' and 'Y' are incompatible. Different number of rows.")
    }
  } else { # if no Y
    if (is.null(pls$data)) 
      stop("Argument 'Y' is missing. No dataset available.")
  }
  if (!is.na(pmatch(method, "bootstrap"))) 
    method <- "bootstrap"
  METHODS <- c("bootstrap", "permutation")
  method <- pmatch(method, METHODS)
  if (is.na(method)) {
    warning("Invalid argument 'method'. Default 'method=bootstrap' is used.")   
    method <- "bootstrap"
  }
  if (is.null(reps) | length(reps)>1) reps<-100
  if (!is.numeric(reps) | floor(reps)<=0) reps<-100

  # =======================================================
  # inputs setting
  # =======================================================  
  IDM <- pls$model$IDM# Inner Design Matrix
  blocks <- pls$model$blocks# cardinality of blocks
  scheme <- pls$model$scheme# inner weighting scheme
  modes <- pls$model$modes# measurement modes
  scaled <- pls$model$scaled# type of scaling
  plsr <- pls$model$plsr# pls-regression
  tol <- pls$model$tol# tolerance criterion
  iter <- pls$model$tol# max num iterations
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
  blocklist <- as.list(1:lvs)
  for (j in 1:lvs)
    blocklist[[j]] <- rep(j,blocks[j])
  blocklist <- unlist(blocklist)
  # apply the selected scaling
  if (scaled) {
    sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
    X <- scale(DM, scale=sd.X)
  } else {
    X <- scale(DM, scale=FALSE)
  }

  # =======================================================
  # Global model estimation
  # =======================================================  
  out.ws <- get_weights(X, IDM, blocks, modes, scheme, tol, iter)
  if (is.null(out.ws)) stop("The pls algorithm is non convergent") 
  cor.XY <- cor(X, X%*%out.ws[[2]])
  w.sig <- rep(NA, lvs)
  for (k in 1:lvs) 
    w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
  if (scaled) {
    Y.lvs <- X %*% out.ws[[2]] %*% diag(w.sig, lvs, lvs)
  } else   
    Y.lvs <- DM %*% out.ws[[2]] %*% diag(w.sig, lvs, lvs)
  dimnames(Y.lvs) = list(rownames(X), lvs.names)
  # Path coefficients
  pathmod <- get_paths(IDM, Y.lvs, plsr)
  innmod <- pathmod[[1]]
  Path.global <- pathmod[[2]]
  R2.global <- pathmod[[3]]
  endo = rowSums(IDM)
  endo[endo != 0] = 1  # vector indicating endogenous LVs
  path.labs = NULL
  efs.labs = NULL
  for (j in 1:lvs)
    for (i in j:lvs)
      if (IDM[i,j]==1) 
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
  if (scaled) {
    sd.Xg1 <- sqrt((ng1-1)/ng1) * apply(DM[group1,], 2, sd)
    X.g1 <- scale(DM[group1,], scale=sd.Xg1)
  } else {
    X.g1 <- scale(DM[group1,], scale=FALSE)
  }
  wgs.g1 <- get_weights(X.g1, IDM, blocks, modes, scheme, tol, iter)
  if (is.null(wgs.g1)) stop("The algorithm is non convergent") 
  cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
  w.sig <- rep(NA, lvs)
  for (k in 1:lvs) 
    w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
  if (scaled) {
    Y1.lvs <- X.g1 %*% wgs.g1[[2]] %*% diag(w.sig, lvs, lvs)
  } else   
    Y1.lvs <- DM[group1,] %*% wgs.g1[[2]] %*% diag(w.sig, lvs, lvs)
  dimnames(Y1.lvs) <- list(rownames(X.g1), lvs.names)
  # Path coefficients 
  pathmod.g1 <- get_paths(IDM, Y1.lvs, plsr)
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
  if (scaled) {
    sd.Xg2 <- sqrt((ng2-1)/ng2) * apply(DM[group2,], 2, sd)
    X.g2 <- scale(DM[group2,], scale=sd.Xg2)
  } else {
    X.g2 <- scale(DM[group2,], scale=FALSE)
  }
  wgs.g2 <- get_weights(X.g2, IDM, blocks, modes, scheme, tol, iter)
  if (is.null(wgs.g2)) stop("The algorithm is non convergent") 
  cor.XY <- cor(X[group2,], X.g2%*%wgs.g2[[2]])
  w.sig <- rep(NA,lvs)
  for (k in 1:lvs) 
    w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
  if (scaled) {
    Y2.lvs <- X.g2 %*% wgs.g2[[2]] %*% diag(w.sig, lvs, lvs)
  } else   
    Y2.lvs <- DM[group2,] %*% wgs.g2[[2]] %*% diag(w.sig, lvs, lvs)
  dimnames(Y2.lvs) <- list(rownames(X.g2), lvs.names)
  # Path coefficients 
  pathmod.g2 <- get_paths(IDM, Y2.lvs, plsr)
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
      if (scaled) {
        sd.Xg1 <- sqrt((ng1-1)/ng1) * apply(DM[samg1,], 2, sd)
        sd.Xg2 <- sqrt((ng2-1)/ng2) * apply(DM[samg2,], 2, sd)
        X.g1 <- scale(DM[samg1,], scale=sd.Xg1)
        X.g2 <- scale(DM[samg2,], scale=sd.Xg2)
      } else {
        X.g1 <- scale(DM[samg1,], scale=FALSE)
        X.g2 <- scale(DM[samg2,], scale=FALSE)
      }
      wgs.g1 <- get_weights(X.g1, IDM, blocks, modes, scheme, tol, iter)
      wgs.g2 <- get_weights(X.g2, IDM, blocks, modes, scheme, tol, iter)
      if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
      if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples") 
      cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
      w.sig <- rep(NA,lvs)
      for (k in 1:lvs) 
        w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
      if (scaled) {
        Y1.lvs <- X.g1 %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs)
      } else   
        Y1.lvs <- DM[samg1,] %*% wgs.g1[[2]] %*% diag(w.sig, lvs, lvs)
      cor.XY <- cor(X.g2, X.g2 %*% wgs.g2[[2]])
      w.sig <- rep(NA, lvs)
      for (k in 1:lvs) 
        w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
      if (scaled) { 
        Y2.lvs <- X.g2 %*% wgs.g2[[2]] %*% diag(w.sig, lvs, lvs)
      } else   
        Y2.lvs <- DM[samg2,] %*% wgs.g2[[2]] %*% diag(w.sig, lvs, lvs)
      pathmod.g1 <- get_paths(IDM, Y1.lvs, plsr)
      paths.g1 <- pathmod.g1[[2]]    
      pathmod.g2 <- get_paths(IDM, Y2.lvs, plsr)
      paths.g2 <- pathmod.g2[[2]]
      BG1[i,] <- as.vector(paths.g1[which(IDM==1)])
      BG2[i,] <- as.vector(paths.g2[which(IDM==1)])
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
      if (scaled) {
        sd.Xg1 <- sqrt((ng1-1)/ng1) * apply(DM[samg1,], 2, sd)
        sd.Xg2 <- sqrt((ng2-1)/ng2) * apply(DM[samg2,], 2, sd)
        X.g1 <- scale(DM[samg1,], scale=sd.Xg1)
        X.g2 <- scale(DM[samg2,], scale=sd.Xg2)
      } else {
        X.g1 <- scale(DM[samg1,], scale=FALSE)
        X.g2 <- scale(DM[samg2,], scale=FALSE)
      }
      wgs.g1 <- get_weights(X.g1, IDM, blocks, modes, scheme, tol, iter)
      wgs.g2 <- get_weights(X.g2, IDM, blocks, modes, scheme, tol, iter)
      if (is.null(wgs.g1)) stop("Non convergence in bootstrap samples") 
      if (is.null(wgs.g2)) stop("Non convergence in bootstrap samples") 
      cor.XY <- cor(X.g1, X.g1%*%wgs.g1[[2]])
      w.sig <- rep(NA, lvs)
      for (k in 1:lvs) 
        w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
      if (scaled) {
        Y1.lvs <- X.g1 %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs)
      } else   
        Y1.lvs <- DM[samg1,] %*% wgs.g1[[2]] %*% diag(w.sig,lvs,lvs)
      cor.XY <- cor(X.g2, X.g2%*%wgs.g2[[2]])
      w.sig <- rep(NA,lvs)
      for (k in 1:lvs) 
        w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
      if (scaled) {
        Y2.lvs <- X.g2 %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs)
      } else   
        Y2.lvs <- DM[samg2,] %*% wgs.g2[[2]] %*% diag(w.sig,lvs,lvs)
      pathmod.g1 <- get_paths(IDM, Y1.lvs, plsr)
      paths.g1 <- pathmod.g1[[2]]    
      pathmod.g2 <- get_paths(IDM, Y2.lvs, plsr)
      paths.g2 <- pathmod.g2[[2]]
      pp1 <- as.vector(paths.g1[which(IDM==1)])
      pp2 <- as.vector(paths.g2[which(IDM==1)])
      dif.perm[i,] <- abs(pp1 - pp2)
    }   
    s.perm <- dif.orig 
    for (i in 1:sum(IDM))         
      s.perm[i] <- length(which(dif.orig[i]<dif.perm[,i])) + 1
    p.val <- (1/(nb+1))*s.perm 
    signi.path <- rep("no",length(p.val))
    signi.path[p.val<0.05] <- "yes"
    res.path <- round(cbind(path.global, path.g1, path.g2, dif.orig, p.val), 4)
    res <- data.frame(res.path, signi.path)
    colnames(res) <- c("global", paste(rep("group",2),levels(g),sep="."), 
                       "diff.abs", "p.value", "sig.05") 
  }

  # =======================================================
  # Results
  # =======================================================  
  met <- switch(method, "1"="bootstrap", "2"="permutation")
  settings <- c(scaled=scaled, scheme=scheme, method=met)
  res = list(test = res, 
             global = innmod, 
             group1 = innmod.g1, 
             group2 = innmod.g2, 
             settings = settings, 
             reps = reps)
  class(res) = "plspm.groups"
  return(res)
}

