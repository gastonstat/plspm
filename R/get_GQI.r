#' @title Group Quality Index
#' @description Calculate Group Quality Index for REBUS-PLS
#' @param pls an object of class \code{"plspm"}
#' @param part vector with unit memberships or categorical variable
#' @param DM data matrix
#' @keywords internal
#' @export
get_GQI <- function(pls, part, DM)
{
  # =========================== arguments ==============================    
  IDM = pls$model$IDM
  blocks = pls$model$blocks
  specs = pls$model$specs
  lvs = pls$model$gens$lvs
  lvs.names = pls$model$gens$lvs_names
  mvs = pls$model$gens$mvs
  blocklist = indexify(blocks)
  endo = rowSums(IDM)
  endo[endo != 0] = 1
  
  # data scaling (standardized data)
  sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
  X <- scale(DM, scale=sd.X)
  clas.part <- part 
  # number of clusters
  nclus <- nlevels(factor(clas.part))
  w.locals <- as.list(1:nclus)# outer.weights
  LV.locals <- as.list(1:nclus)# std latent variables 
  loads.locals <- as.list(1:nclus)# loadings
  path.locals <- as.list(1:nclus)# path coefficients
  R2.locals <- as.list(1:nclus)# R2
  comu.locals <- as.list(1:nclus)# mvs communalities
  outres.locals <- as.list(1:nclus)# communality residuals
  innres.locals <- as.list(1:nclus)# structural residuals
  out.term <- as.list(1:nclus)# outer term for GQI
  inn.term <- as.list(1:nclus)# inner term for GQI
  gqi.locals <- rep(0, nclus)# pseudo-gqi for each class
  
  # define MV matrix for each initial class
  split.DM <- as.list(1:nclus)
  split.X <- as.list(1:nclus)
  for (k in 1:nclus)
    split.DM[[k]] <- DM[clas.part==k,]            
  
  # local models computation
  for (k in 1:nclus)
  {   
    nk <- nrow(split.DM[[k]])
    mean.k <- apply(split.DM[[k]],2,mean)# local mean
    sd.k <- sqrt((nk-1)/nk) * apply(split.DM[[k]],2,sd)# local std.dev
    # spliting data matrix for each class
    split.X[[k]] <- scale(split.DM[[k]], center=mean.k, scale=sd.k)
    # calculating outer weights for each class
    out.ws  <- get_weights(split.X[[k]], IDM, blocks, specs)
    w.locals[[k]] <- out.ws$W
    # calculating LV scores for each class
    LV.locals[[k]] <- split.X[[k]] %*% out.ws$W
    # calculating path coefficients for each class
    pathmod <- get_paths(IDM, LV.locals[[k]])
    path.locals[[k]] <- pathmod[[2]]
    R2.locals[[k]] <- pathmod[[3]][endo==1]
    # calculating loadings and communalities for each class
    xloads = cor(split.X[[k]], LV.locals[[k]])
    loads.locals[[k]] = rowSums(xloads * out.ws$ODM)
    comu.locals[[k]] = loads.locals[[k]]^2
    # computation of communality residuals (squared)
    out.res <- split.X[[k]]
    for (j in 1:lvs)
    {
      X.hat <- LV.locals[[k]][,j] %*% t(loads.locals[[k]][blocklist==j])
      # outer residuals
      out.res[,blocklist==j] <- (split.X[[k]][,blocklist==j] - X.hat)^2
    }
    outres.locals[[k]] <- out.res
    # computation of inner residuals (squared)
    if (sum(endo)!=1)
      Y.hat <- LV.locals[[k]] %*% t(path.locals[[k]][endo==1,])   
    if (sum(endo)==1)
      Y.hat <- LV.locals[[k]] %*% path.locals[[k]][endo==1,]        
    innres.locals[[k]] <- (LV.locals[[k]][,endo==1] - Y.hat)^2
    # outer and inner terms of GQI formula
    aux_out_term = 1 - (colSums(outres.locals[[k]]) / colSums(split.X[[k]]^2))
    out.term[[k]] <- mean(aux_out_term)
    if (sum(endo) == 1) {
      inn_numer = colSums(innres.locals[[k]])
      inn_denom = sum(LV.locals[[k]][,endo==1]^2)
      inn.term[[k]] = mean(1 - (inn_numer / inn_denom))
    }
    if (sum(endo) != 1) {
      inn_numer = colSums(innres.locals[[k]])
      inn_denom = colSums(LV.locals[[k]][,endo==1]^2)
      inn.term[[k]] = mean(1 - (inn_numer / inn_denom))
    }
    gqi.locals[k] <- out.term[[k]] * inn.term[[k]]
  }
  # proportion of units in each class
  unit.prop = unlist(lapply(split.X, nrow)) / nrow(DM)
  GQI = sqrt(sum(gqi.locals * unit.prop))
  return(GQI)
}
