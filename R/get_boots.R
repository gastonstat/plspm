#' @title Performs bootstrap validation in \code{plspm}
#' 
#' @description
#' Internal function. \code{get_boots} is called by \code{plspm}.
#' 
#' @param DM Data Matrix
#' @param IDM Inner Design Matrix
#' @param blocks vector with number of lvs per block
#' @param modes vector of modes
#' @param scheme type of inner path scheme
#' @param scaled logical to scale the data
#' @param br number of bootstrap resamples
#' @param plsr logical to get path coefficients by pls regression
#' @param tol threshold for convergence
#' @param iter maximum number of iterations
#' @keywords internal
#' @export
get_boots <-
  function(DM, IDM, blocks, modes, scheme, scaled, br, plsr, tol, iter)
  {
    # =======================================================
    # inputs setting
    # =======================================================  
    lvs <- nrow(IDM)
    lvs.names <- rownames(IDM)
    mvs <- ncol(DM)
    mvs.names <- colnames(DM)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
      blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1    
    bootnum <- br
    # scaling data
    if (scaled) {
      sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
      X <- scale(DM, scale=sd.X)
    } else {
      X <- scale(DM, scale=FALSE)
    }
    colnames(X) <- mvs.names

    # =======================================================
    # computation of the original plspm model
    # =======================================================  
    out.ws <- get_weights(X, IDM, blocks, modes, scheme, tol, iter)
    wgs.orig <- out.ws[[1]]
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
      w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Y.lvs <- X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
    pathmod <- get_paths(IDM, Y.lvs, plsr)
    Path <- pathmod[[2]]
    path.orig <- as.vector(Path[which(IDM==1)])
    r2.orig <- pathmod[[3]][which(endo==1)]
    Path.efs <- get_effects(Path)
    loadcomu <- get_loads(X, Y.lvs, blocks)    
    loads.orig <- loadcomu[[1]]

    # =======================================================
    # Bootstrap Validation
    # =======================================================  
    path.labs <- NULL
    efs.labs <- NULL
    for (j in 1:lvs)
      for (i in j:lvs)
        if (IDM[i,j]==1) 
          path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))    
    WEIGS <- matrix(NA, bootnum, mvs)
    LOADS <- matrix(NA, bootnum, mvs)
    PATHS <- matrix(NA, bootnum, sum(IDM))
    TOEFS <- matrix(NA, bootnum, nrow(Path.efs))
    RSQRS <- matrix(NA, bootnum, sum(endo))
    i <- 1
    while (i <= bootnum)
    {
      boot.obs <- sample.int(nrow(X), size=nrow(X), replace=TRUE)
      DM.boot <- DM[boot.obs,]
      # scaling boot sample
      if (scaled) {
        sd.XB <- sqrt((nrow(DM.boot)-1)/nrow(DM.boot)) * apply(DM.boot, 2, sd)
        X.boot <- scale(DM.boot, scale=sd.XB)
      } else {
        X.boot <- scale(DM.boot, scale=FALSE)
      }
      colnames(X.boot) <- mvs.names
      # calculating boot model parameters 
      w.boot <- get_weights(X.boot, IDM, blocks, modes, scheme, tol, iter)
      if (is.null(w.boot)) {
        i <- i - 1
        next
      }
      WEIGS[i,] <- w.boot[[1]]
      Y.boot <- X.boot %*% w.boot[[2]]
      pathmod <- get_paths(IDM, Y.boot, plsr)
      P.boot <- pathmod[[2]]
      Toef.boot <- get_effects(P.boot)
      PATHS[i,] <- as.vector(P.boot[which(IDM==1)])
      TOEFS[i,] <- Toef.boot[,4]
      RSQRS[i,] <- pathmod[[3]][which(endo==1)]
      l.boot <- get_loads(X.boot, Y.boot, blocks)    
      LOADS[i,] <- l.boot[[1]]
      i <- i + 1
    }

    # =======================================================
    # Bootstrap results
    # =======================================================  
    # Outer weights
    colnames(WEIGS) <- mvs.names
    WB <- data.frame(Original = wgs.orig, 
                     Mean.Boot = apply(WEIGS, 2, mean), 
                     Std.Error = apply(WEIGS, 2, sd), 
                     perc.025 = apply(WEIGS, 2, function(x) quantile(x, 0.025)),
                     perc.975 = apply(WEIGS, 2, function(x) quantile(x, 0.975)))
    # Loadings
    colnames(LOADS) <- mvs.names
    LB <- data.frame(Original = loads.orig, 
                     Mean.Boot = apply(LOADS, 2, mean),
                     Std.Error = apply(LOADS, 2, sd), 
                     perc.025 = apply(LOADS, 2, function(x) quantile(x, 0.025)),
                     perc.975 = apply(LOADS, 2, function(x) quantile(x, 0.975)))
    # Path coefficients
    colnames(PATHS) <- path.labs
    PB <- data.frame(Original = path.orig, 
                     Mean.Boot = apply(PATHS, 2, mean),
                     Std.Error = apply(PATHS, 2, sd), 
                     perc.025 = apply(PATHS, 2, function(x) quantile(x, 0.025)),
                     perc.975 = apply(PATHS, 2, function(x) quantile(x, 0.975)))
    # Total effects
    colnames(TOEFS) <- Path.efs[, 1]
    TE <- data.frame(Original = Path.efs[, 4], 
                     Mean.Boot = apply(TOEFS, 2, mean), 
                     Std.Error = apply(TOEFS, 2, sd),
                     perc.025 = apply(TOEFS, 2, function(x) quantile(x, 0.025)), 
                     perc.975 = apply(TOEFS, 2, function(x) quantile(x, 0.975)))
    # R-squared
    colnames(RSQRS) <- lvs.names[endo == 1]
    RB <- data.frame(Original = r2.orig, 
                     Mean.Boot = apply(RSQRS, 2, mean),
                     Std.Error = apply(RSQRS, 2, sd), 
                     perc.025 = apply(RSQRS, 2, function(x) quantile(x, 0.025)),
                     perc.975 = apply(RSQRS, 2, function(x) quantile(x, 0.975)))
    # Bootstrap Results
    res.boot <- list(weights = WB, 
                     loadings = LB, 
                     paths = PB, 
                     rsq = RB, 
                     total.efs = TE)
    return(res.boot)
  }
