#' @title Local groups comparison test
#' 
#' @description Calculate Group Quality Index for REBUS-PLS
#' Internal function.
#' 
#' @param X data matrix
#' @param pls an object of class \code{"plspm"}
#' @param g factor with 2 levels indicating the groups to be compared
#' @return list with path coefficients, loadings, and gof
#' @keywords internal
#' @export
get_locals_test <- function(X, pls, g)
{
  # ========================== INPUTS SETTING ==========================
  IDM = pls$model$IDM
  blocks = pls$model$blocks
  specs = pls$model$specs
  scaled = pls$model$specs$scaled
  
  lvs = pls$model$gens$lvs
  lvs.names = pls$model$gens$lvs_names
  mvs = pls$model$gens$mvs
  blocklist = indexify(blocks)
  reps <- 100
  path.labs <- NULL
  efs.labs <- NULL
  for (j in 1:lvs)
    for (i in j:lvs)
      if (IDM[i,j]==1) 
        path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
  
  # ====================== Group1 model estimation =====================
  g1.lab <- levels(g)[1]
  group1 <- which(g==levels(g)[1])
  ng1 <- length(group1)
  # apply the selected scaling
  X.g1 = get_data_scaled(X[group1,], scaled)
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
  # calculating loadings
  xloads = cor(X.g1, Y1.lvs)
  load.g1 = rowSums(xloads * wgs.g1$ODM)
  # gof
  gof.g1 = get_gof(load.g1^2, R2.g1, blocks, IDM)
  
  # ====================== Group2 model estimation =====================
  g2.lab <- levels(g)[2]
  group2 <- which(g==levels(g)[2])
  ng2 <- length(group2)
  # apply the selected scaling
  X.g2 = get_data_scaled(X[group2,], scaled)
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
  # calculating loadings and communalities for each class
  xloads = cor(X.g2, Y2.lvs)
  load.g2 = rowSums(xloads * wgs.g1$ODM)
  # gof
  gof.g2 = get_gof(load.g2^2, R2.g2, blocks, IDM)
  
  # ====================== Group Comparison =====================
  difpath.orig = abs(path.g1 - path.g2)
  difload.orig = abs(load.g1 - load.g2)
  difgof.orig = abs(gof.g1 - gof.g2)
  group1 = which(g==levels(g)[1])
  group2 = which(g==levels(g)[2])
  ng1 = length(group1)
  ng2 = length(group2)
  difpath.perm = matrix(0, reps, sum(IDM))
  difload.perm = matrix(0, reps, mvs)
  difgof.perm = rep(0, reps)
  for (i in 1:reps)# multigroup permutation
  {
    permu <- sample(1:(ng1+ng2), ng1+ng2)
    samg1 <- permu[1:ng1]
    samg2 <- permu[(ng1+1):(ng1+ng2)]
    
    # apply the selected scaling
    X.g1 = get_data_scaled(X[samg1,], scaled)
    wgs.g1 = get_weights(X.g1, IDM, blocks, specs)
    ok_weights = test_null_weights(wgs.g1, specs)  
    Y1.lvs = get_scores(X.g1, wgs.g1$W)
    
    # apply the selected scaling
    X.g2 = get_data_scaled(X[samg1,], scaled)
    wgs.g2 = get_weights(X.g2, IDM, blocks, specs)
    ok_weights = test_null_weights(wgs.g2, specs)  
    Y2.lvs = get_scores(X.g2, wgs.g2$W)
    
    # Path coefficients 
    pathmod.g1 <- get_paths(IDM, Y1.lvs)
    paths.g1 <- pathmod.g1[[2]]
    R2.g1 <- pathmod.g1[[3]]    
    pathmod.g2 <- get_paths(IDM, Y2.lvs)
    paths.g2 <- pathmod.g2[[2]]
    R2.g2 <- pathmod.g2[[3]]    

    # calculating loadings
    xloads = cor(X.g1, Y1.lvs)
    load.g1 = rowSums(xloads * wgs.g1$ODM)
    xloads = cor(X.g2, Y2.lvs)
    load.g2 = rowSums(xloads * wgs.g2$ODM)

    # gof
    gof.g1 = get_gof(load.g1^2, R2.g1, blocks, IDM)
    gof.g2 = get_gof(load.g2^2, R2.g2, blocks, IDM)
    
    # difference between groups
    pp1 <- as.vector(paths.g1[which(IDM==1)])
    pp2 <- as.vector(paths.g2[which(IDM==1)])
    difpath.perm[i,] <- abs(pp1 - pp2)
    difload.perm[i,] <- abs(load.g1 - load.g2)
    difgof.perm[i] <- abs(gof.g1 - gof.g2)
  }   
  # p-value for path coefficients
  path.perm <- difpath.orig 
  for (j in 1:sum(IDM))         
    path.perm[j] <- length(which(difpath.orig[j]<difpath.perm[,j])) + 1
  path.val <- (1/(reps+1)) * path.perm 
  signi.path <- rep("no", length(path.val))
  signi.path[path.val < 0.05] <- "yes"
  res.path <- round(cbind(path.g1, path.g2, difpath.orig, path.val), 4)
  res1 <- data.frame(res.path, signi.path)
  colnames(res1) <- c(paste(rep("Class",2), levels(g),sep="."), 
                      "diff.abs", "p.value", "sig.05")  
  # p-values for loadings
  load.perm <- difload.orig 
  for (j in 1:mvs)
    load.perm[j] <- length(which(difload.orig[j]<difload.perm[,j])) + 1
  load.val <- (1/(reps+1)) * load.perm 
  signi.load <- rep("no", length(load.val))
  signi.load[load.val < 0.05] <- "yes"
  res.load <- round(cbind(load.g1, load.g2, difload.orig, load.val), 4)
  res2 <- data.frame(res.load, signi.load)
  colnames(res2) <- c(paste(rep("Class",2), levels(g),sep="."), 
                      "diff.abs", "p.value", "sig.05")  
  # p-values for gof
  gof.perm <- length(which(difgof.orig<difgof.perm)) + 1
  gof.val <- (1/(reps+1)) * gof.perm 
  signi.gof <- rep("no", length(gof.val))
  signi.gof[gof.val < 0.05] <- "yes"
  res3 <- data.frame(round(gof.g1,4), round(gof.g2,4), 
                     round(difgof.orig,4), round(gof.val,4), signi.gof)
  names(res3) <- c(paste(rep("Class",2), levels(g),sep="."), 
                   "diff.abs", "p.value", "sig.05")  
  # list with results 
  resul <- list(paths=res1, loadings=res2, gof=res3)
  return(resul)
}
