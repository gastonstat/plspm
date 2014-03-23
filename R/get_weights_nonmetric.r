#' @title Outer Weights Non-Metric Data
#' 
#' @details
#' Internal function. \code{get_weights} is called by \code{plspm}
#' 
#' @note
#' Calculate outer weights for non-metric data (under Lohmoller's algorithm)
#' 
#' @param X scaled data
#' @param path_matrix matrix with path connections
#' @param blocks list with variables in each block
#' @param specs list with algorithm specifications
#' @return list of outer weights, ODM, scores, QQ, iter
#' @keywords internal
#' @template internals
#' @export
get_weights_nonmetric <-
function(X, path_matrix, blocks, specs)
{
  lvs = nrow(path_matrix)
  mvs = ncol(X)
  num_obs = nrow(X)
  correction = sqrt(nrow(X) / (nrow(X)-1))
  blockinds = indexify(blocks)
  block_sizes = lengths(blocks)
  PLScomp = specs$plscomp
  start_end = from_to(block_sizes)
  
  # create dummy matrices for categorical manifest variables
  dummies = get_dummies(X, specs)
  
  ## transforming X in a list of blocks
  Xblocks = vector("list", lvs)
  start_end = from_to(blocks)
  from = start_end$from
  to = start_end$to
  for (q in 1:lvs) {
    if (from[q] == to[q]) {
      Xblocks[[q]] = as.matrix(X[,from[q]:to[q]])
    } else {
      Xblocks[[q]] = X[,from[q]:to[q]]      
    }
  }
  
  # list for quantification of blocks' variables
  QQ = Xblocks
  # missing data flags
  missing_data = sapply(Xblocks, is_missing)
  # initialize list with availibility indicators (to handle NAs)
  X_avail = vector("list", lvs)
  # outer design matrix 'ODM' and matrix of outer weights 'W'
  ODM = list_to_dummy(blocks)
  
  # =======================================================================
  # initialization
  # =======================================================================
  # outer weights (normalized)
  w_ones = list_ones(block_sizes)
  w = lapply(w_ones, normalize)
  # LV scores
  Y = matrix(0, num_obs, lvs)
  for (q in 1:lvs) {
    if (missing_data[q]) {
      # binary matrix (1=available data, 0=NA)
      X_avail[[q]] = 1 - is.na(Xblocks[[q]])
      for (i in 1:nrow(X)) {
        aux_numerator = sum(QQ[[q]][i,]*w[[q]], na.rm = TRUE)
        aux_denom = sum(w[[q]][which(is.na(QQ[[q]][i,]*w[[q]]) == FALSE)]^2)
        Y[i,q] <- aux_numerator / aux_denom
      }
    } else {
      Y[,q] = QQ[[q]] %*% w[[q]]        
    }
  }
  
  # =======================================================================
  # iterative cycle
  # =======================================================================
  # matrix of inner weights
  E = matrix(0, lvs, lvs)
  link = t(path_matrix) + path_matrix
  z_temp = matrix(0, num_obs, 1)
  iter = 0
  repeat 
  {
    iter = iter + 1
    #    y_old = as.numeric(Y)
    Y_old = Y
    
    # =============================================================
    # updating inner weights
    # =============================================================
    E <- switch(specs$scheme, 
                "centroid" = sign(cor(Y) * link),
                "factorial" = cov(Y) * link,
                "path" = get_path_scheme(path_matrix, Y))
    # internal estimation of LVs 'Z'
    Z = Y %*% E
    
    # for each block
    for (q in 1:lvs) 
    {
      # standardize inner estimates if PLScore mode
      # if (specs$modes[q] != "PLSCOW") {
      #### Giorgio's suggestion: do not standardize inner estimates:
       #if (specs$modes[q] != "PLSCOW" & specs$modes[q] != "NEWA") {
       #  Z[,q] <- scale(Z[,q]) * correction
       #}
      # =============================================================
      # Quantification of the MVs in block ["QQ"]
      # =============================================================
      # for each MV in block 'q'
      if (specs$modes[q] == "B" && block_sizes[q] > 1) {
        Beta <- summary(lm(Z[,q]~QQ[[q]]))$coef[-1,1]
        X.star <- matrix(,num_obs,block_sizes[q])
        for (p in 1L:block_sizes[q]) {
          X.star[,p] <- (1/Beta[p])*(Z[,q] - (QQ[[q]][,-p,drop=FALSE]%*%Beta[-p]))
          if (specs$scaling[[q]][p] == "nom") {
            # extract corresponding dummy matrix
            #          aux_dummy = dummies[[blocks[[q]][p]]]
            which_dummy = (start_end$from[q]:start_end$to[q])[p]
            aux_dummy = dummies[[which_dummy]]
            # apply scaling
            QQ[[q]][,p] = get_nom_scale(X.star[,p], Xblocks[[q]][,p], aux_dummy)
            QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
          }
          if (specs$scaling[[q]][p] == "ord") {
            # extract corresponding dummy matrix
            #          aux_dummy = dummies[[blocks[[q]][p]]]
            which_dummy = (start_end$from[q]:start_end$to[q])[p]
            aux_dummy = dummies[[which_dummy]]
            # apply scaling
            QQ[[q]][,p] = get_ord_scale(X.star[,p], Xblocks[[q]][,p], aux_dummy)
            QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
          }                   
          if (specs$scaling[[q]][p] == "num") {
            QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
          }
        }
      }
      else {
        for (p in 1L:block_sizes[q]) {
          if (specs$scaling[[q]][p] == "nom") {
            # extract corresponding dummy matrix
            #          aux_dummy = dummies[[blocks[[q]][p]]]
            which_dummy = (start_end$from[q]:start_end$to[q])[p]
            aux_dummy = dummies[[which_dummy]]
            # apply scaling
            QQ[[q]][,p] = get_nom_scale(Z[,q], Xblocks[[q]][,p], aux_dummy)
            QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
          }
          if (specs$scaling[[q]][p] == "ord") {
            # extract corresponding dummy matrix
            #          aux_dummy = dummies[[blocks[[q]][p]]]
            which_dummy = (start_end$from[q]:start_end$to[q])[p]
            aux_dummy = dummies[[which_dummy]]
            # apply scaling
            QQ[[q]][,p] = get_ord_scale(Z[,q], Xblocks[[q]][,p], aux_dummy)
            QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
          }                   
          if (specs$scaling[[q]][p] == "num") {
            QQ[[q]][,p] = get_num_scale(QQ[[q]][,p])
          }
          ### DO WE REALLY NEED THIS LINE:
          #if (specs$scaling[[q]][p] == "raw") {
          #  QQ[[q]][,p] = QQ[[q]][,p]
          #}
        }
      }
      
      # =============================================================
      # updating outer weights "w" and outer estimates "Y"
      # =============================================================
      
      # Mode A (="PLScore1comp") ====================================
      if (specs$modes[q] == "A") {
        if (missing_data[q]) {
          # compute w[[q]][l] as the regr. coeff. of QQ[[q]][,l] on Z[,q] 
          # considering only the lines where QQ[[q]][i,l] exist
          w[[q]] = colSums(QQ[[q]]*Z[,q], na.rm = TRUE)
          w[[q]] = w[[q]] / colSums((X_avail[[q]]*Z[,q])^2)
          # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
          # considering only the columns where QQ[[q]][i,l] exist
          Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm = TRUE)
          Y[,q] = Y[,q] / colSums((t(X_avail[[q]])*w[[q]])^2)
          # normalize Y[,q] to unitary variance
          Y[,q] = scale(Y[,q]) * correction
        } 
        else {# complete data in block q
          w[[q]] = (t(QQ[[q]]) %*% Z[,q]) / sum(Z[,q]^2)
          Y[,q] = QQ[[q]] %*% w[[q]]
          Y[,q] = scale(Y[,q]) * correction
        }
      }
      
      #  Mode New A (= "PLScow1comp") ================================
      if (specs$modes[q] == "NEWA") {
        if (missing_data[q]) {
          # compute w[[q]][l] as the regr. coeff. of QQ[[q]][,l] on Z[,q]
          # considering just the lines where QQ[[q]][i,l] exist
          w[[q]] = colSums(QQ[[q]]*Z[,q], na.rm = TRUE)
          w[[q]] = w[[q]]/colSums((X_avail[[q]]*Z[,q])^2)
          # normalize w[[q]] to unitary norm
          w[[q]] = w[[q]] / sqrt(sum(w[[q]]^2))
          # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
          # considering only the columns where QQ[[q]][i,l] exist
          Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm=TRUE)
          Y[,q] = Y[,q] / colSums((t(X_avail[[q]])*w[[q]])^2)
        } 
        else {# complete data in block q
          w[[q]] = (t(QQ[[q]]) %*% Z[,q]) / sum(Z[,q]^2)
          w[[q]] = w[[q]] / sqrt(sum(w[[q]]^2))
          Y[,q] = QQ[[q]] %*% w[[q]]
        }
      }
      
      # Mode B (NAs were not allowed.. so far. Now we can use PLSR) ====
      if (specs$modes[q] == "B") {
        if (missing_data[q]) {# use full component PLS-R 
          w[[q]] = get_PLSR_NA(Y = Z[,q], X = QQ[[q]], ncomp = block_sizes[q])$B
          # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
          # considering only the columns where QQ[[q]][i,l] exist
          Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm=TRUE)
          Y[,q] = Y[,q]/colSums((t(X_avail[[q]])*w[[q]])^2)
          # normalize Y[,q] to unitary variance
          Y[,q] = scale(Y[,q]) * correction		
        }
        else {# complete data in block q
          w[[q]] = solve.qr(qr(QQ[[q]]), Z[,q])
          #w[[q]] = solve(t(QQ[[q]]) %*% QQ[[q]]) %*% t(QQ[[q]]) %*% Z[,q]
          Y[,q] = QQ[[q]] %*% w[[q]]
          Y[,q] = scale(Y[,q]) * correction
        }	
      }
      
      # Mode PLScore ===================================================
      if (specs$modes[q] == "PLSCORE") {
        if (missing_data[q]) {
          w[[q]] = get_PLSR_NA(Y = Z[,q], X = QQ[[q]], ncomp = PLScomp[q])$B
          # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
          # considering only the columns where QQ[[q]][i,l] exist
          Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm=TRUE)
          Y[,q] = Y[,q]/colSums((t(X_avail[[q]])*w[[q]])^2)
          # normalize Y[,q] to unitary variance
          Y[,q] = scale(Y[,q]) * correction		
        }
        else {# complete data in block q
          w[[q]] = get_PLSR(Y = Z[,q], X = QQ[[q]], ncomp = PLScomp[q])$B
          Y[,q] = QQ[[q]] %*% w[[q]]
          Y[,q] = scale(Y[,q]) * correction
        }	
      }
      
      # Mode PLScow =====================================================
      if (specs$modes[q] == "PLSCOW") {
        if (missing_data[q]) {
          w[[q]] = get_PLSR_NA(Y = Z[,q], X = QQ[[q]], ncomp = PLScomp[q])$B
          # normalize w[[q]] to unitary norm
          w[[q]] = w[[q]] / sqrt(sum(w[[q]]^2))
          # compute Y[i,q] as the regr. coeff. of QQ[[q]][i,] on w[[q]] 
          # considering only the columns where QQ[[q]][i,l] exist
          Y[,q] = colSums(t(QQ[[q]])*w[[q]], na.rm=TRUE)
          Y[,q] = Y[,q]/colSums((t(X_avail[[q]])*w[[q]])^2)
        }
        else {# complete data in block q
          w[[q]] = get_PLSR(Y = Z[,q], X = QQ[[q]], ncomp = PLScomp[q])$B
          w[[q]] = w[[q]] / sqrt(sum(w[[q]]^2))
          Y[,q] = QQ[[q]] %*% w[[q]]
        }	
      }
      
    }
    # check convergence
    convergence <- sum((abs(Y_old) - abs(Y))^2)
    # Y_old: keep it as a matrix
    #    convergence <- sum((abs(y_old) - abs(as.numeric(Y)))^2)
    if (convergence < specs$tol | iter > specs$maxiter) 
      break
  } # end repeat
  
  # preparing results
  if (iter == specs$maxiter) {
    results = NULL
  } else {
    W = list_to_matrix(lapply(w, as.numeric))
    # open new lines
    QQ = do.call("cbind", QQ)
    W = W %*% diag(1/(apply(QQ %*% W, 2, sd, na.rm=TRUE)/correction), lvs, lvs)
    # end new lines
    w = rowSums(W)
    dimnames(W) = list(colnames(X), rownames(path_matrix))
    dimnames(Y) = list(rownames(X), rownames(path_matrix))
    results = list(w = w, W = W, Y = Y, QQ = QQ, ODM = ODM, iter = iter)
  }
  # output
  results  
}
