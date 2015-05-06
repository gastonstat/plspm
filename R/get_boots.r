#' @title Performs bootstrap validation in \code{plspm}
#' 
#' @details
#' Internal function. \code{get_boots} is called by \code{plspm}.
#' 
#' @param DM Data Matrix
#' @param path_matrix Inner Design Matrix
#' @param blocks list of vectors with column indices
#' @param specs list with algorithm specifications
#' @param br number of bootstrap resamples
#' @keywords internal
#' @template internals
#' @export
get_boots <-
function(DM, path_matrix, blocks, specs, br)
{
  # =======================================================
  # inputs setting
  # =======================================================  
  lvs = nrow(path_matrix)
  lvs.names = rownames(path_matrix)
  mvs = ncol(DM)
  mvs.names = colnames(DM)
  blocklist = indexify(blocks)
  endo = sign(rowSums(path_matrix))
  bootnum = br
  # apply corresponding treatment (centering, reducing, ranking)
  X = get_treated_data(DM, specs)
  
  # =======================================================
  # computation of the original plspm model
  # =======================================================  
  metric = get_metric(specs$scaling)
  if (metric) {
    # object 'weights' contains outer w's, W, ODM, iter
    out.ws = get_weights(X, path_matrix, blocks, specs)
    ok_weights = test_null_weights(out.ws, specs)
    wgs.orig = out.ws$w
    Y.lvs = get_scores(X, out.ws$W)
  } else {
    # object 'weights' contains outer w's, W, Y, QQ, ODM, iter
    out.ws = get_weights_nonmetric(X, path_matrix, blocks, specs)
    ok_weights = test_null_weights(out.ws, specs)
    wgs.orig = out.ws$w
    Y.lvs = out.ws$Y
    X = out.ws$QQ  # quantified MVs
  }
  
  pathmod <- get_paths(path_matrix, Y.lvs)
  Path <- pathmod[[2]]
  path.orig <- as.vector(Path[path_matrix==1])
  r2.orig <- pathmod[[3]][endo==1]
  Path.efs <- get_effects(Path)
  xloads = cor(X, Y.lvs)
  loads.orig = rowSums(xloads * out.ws$ODM)
  
  # =======================================================
  # Bootstrap Validation
  # =======================================================  
  path.labs <- NULL
  efs.labs <- NULL
  for (j in 1:lvs)
    for (i in j:lvs)
      if (path_matrix[i,j]==1) 
        path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i]))
  WEIGS <- matrix(NA, bootnum, mvs)
  LOADS <- matrix(NA, bootnum, mvs)
  PATHS <- matrix(NA, bootnum, sum(path_matrix))
  TOEFS <- matrix(NA, bootnum, nrow(Path.efs))
  RSQRS <- matrix(NA, bootnum, sum(endo))
  i <- 1
  while (i <= bootnum)
  {
    boot.obs <- sample.int(nrow(X), size=nrow(X), replace=TRUE)
    DM.boot <- DM[boot.obs,]
    # apply corresponding treatment (centering, reducing, ranking)
    X.boot = get_treated_data(DM.boot, specs)        
    # calculating boot model parameters 
    if (metric) {
      # object 'weights' contains outer w's, W, ODM, iter
      w.boot = get_weights(X.boot, path_matrix, blocks, specs)
      if (is.null(w.boot)) {
        i <- i - 1
        next
      }
      Y.boot = get_scores(X.boot, w.boot$W)
    } else {
      # object 'weights' contains outer w's, W, Y, QQ, ODM, iter
      w.boot = get_weights_nonmetric(X.boot, path_matrix, blocks, specs)
      if (is.null(w.boot)) {
        i <- i - 1
        next
      }
      Y.boot = w.boot$Y
      X.boot = w.boot$QQ  # quantified MVs
      # X.boot = do.call(cbind, w.boot$QQ)  # quantified MVs
    }
    WEIGS[i,] <- w.boot$w
    pathmod <- get_paths(path_matrix, Y.boot)
    P.boot <- pathmod[[2]]
    Toef.boot <- get_effects(P.boot)
    PATHS[i,] <- as.vector(P.boot[path_matrix==1])
    TOEFS[i,] <- Toef.boot[,4]
    RSQRS[i,] <- pathmod[[3]][endo==1]
    xloads = cor(X.boot, Y.boot)
    LOADS[i,] = rowSums(xloads * w.boot$ODM)
    i <- i + 1
  }
  
  # =======================================================
  # Bootstrap results
  # ======================================================= 
  # Outer weights
  WB = get_boot_stats(WEIGS, wgs.orig)
  #rownames(WB) = mvs.names
  rownames(WB) <- paste(rep(lvs.names, sapply(blocks, length)), 
                        mvs.names, sep='-')
  # Loadings
  LB = get_boot_stats(LOADS, loads.orig)
  #rownames(LB) = mvs.names
  rownames(LB) <- paste(rep(lvs.names, sapply(blocks, length)), 
                       mvs.names, sep='-')
  # Path coefficients
  colnames(PATHS) = path.labs
  PB = get_boot_stats(PATHS, path.orig)
  # R-squared
  colnames(RSQRS) = lvs.names[endo == 1]
  RB = get_boot_stats(RSQRS, r2.orig)
  # Total effects
  colnames(TOEFS) = Path.efs[, 1]
  TE = get_boot_stats(TOEFS, Path.efs[,4]) 
  
  # Bootstrap Results
  list(weights = WB, 
       loadings = LB, 
       paths = PB, 
       rsq = RB, 
       total.efs = TE)
}


#' @title Get data frame with bootstrap statistics
#' 
#' @details
#' Internal function. \code{get_boot_stats} is called by \code{get_boots}.
#' 
#' @param MATRIX Matrix with bootstrapped results
#' @param original vector with original values
#' @keywords internal
#' @template internals
#' @export
get_boot_stats <- function(MATRIX, original) {
  data.frame(Original = original,
             Mean.Boot = apply(MATRIX, 2, mean), 
             Std.Error = apply(MATRIX, 2, sd), 
             perc.025 = apply(MATRIX, 2, function(x) quantile(x, 0.025)),
             perc.975 = apply(MATRIX, 2, function(x) quantile(x, 0.975)))
}
