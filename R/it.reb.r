#' @title Iterative steps of Response-Based Unit Segmentation (REBUS)
#' 
#' @description
#' REBUS-PLS is an iterative algorithm for performing response based 
#' clustering in a PLS-PM framework. \code{it.reb} allows to perform the
#' iterative steps of the REBUS-PLS Algorithm. 
#' It provides summarized results for final local models and the final 
#' partition of the units. Before running this function, it is necessary 
#' to run the \code{\link{res.clus}} function to choose the number of 
#' classes to take into account.
#' 
#' @param pls an object of class \code{"plspm"}
#' @param hclus.res object of class \code{"res.clus"} returned 
#' by \code{\link{res.clus}}
#' @param nk integer larger than 1 indicating the number of classes.
#' This value should be defined according to the dendrogram obtained
#' by performing \code{\link{res.clus}}.
#' @param Y optional data matrix used when \code{pls$data} is \code{NULL}
#' @param stop.crit Number indicating the stop criterion for the iterative 
#' algorithm. It is suggested to use the threshold of less than 0.05\% 
#' of units changing class from one iteration to the other as stopping rule.
#' @param iter.max integer indicating the maximum number of iterations
#' @return an object of class \code{"rebus"}
#' @return \item{loadings}{Matrix of standardized loadings (i.e. correlations 
#' with LVs.) for each local model}
#' @return \item{path.coefs}{Matrix of path coefficients for each local model}
#' @return \item{quality}{Matrix containing the average communalities, the
#' average redundancies, the R2 values, and the GoF index for each local model}
#' @return \item{segments}{Vector defining the class membership of each unit} 
#' @return \item{origdata.clas}{ The numeric matrix with original data and with 
#' a new column defining class membership of each unit}
#' @author Laura Trinchera, Gaston Sanchez
#' 
#' @references Esposito Vinzi, V., Trinchera, L., Squillacciotti, S., 
#' and Tenenhaus, M. (2008) REBUS-PLS: A Response-Based Procedure for detecting 
#' Unit Segments in PLS Path Modeling. \emph{Applied Stochastic Models in 
#' Business and Industry (ASMBI)}, \bold{24}, pp. 439-458.
#' 
#' Trinchera, L. (2007) Unobserved Heterogeneity in Structural Equation Models: 
#' a new approach to latent class detection in PLS Path Modeling. 
#' \emph{Ph.D. Thesis}, University of Naples "Federico II", Naples, Italy.
#' 
#' \url{http://www.fedoa.unina.it/2702/1/Trinchera_Statistica.pdf}
#' @seealso \code{\link{plspm}}, \code{\link{rebus.pls}}, 
#' \code{\link{res.clus}}
#' @export
#' @examples
#' \dontrun{
#' ## Example of REBUS PLS with simulated data
#' 
#' # load simdata
#' data("simdata", package='plspm')
#' 
#' # Calculate global plspm
#' sim_inner = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
#' dimnames(sim_inner) = list(c("Price", "Quality", "Satisfaction"),
#'                            c("Price", "Quality", "Satisfaction"))
#' sim_outer = list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13)) 
#' sim_mod = c("A", "A", "A")  # reflective indicators
#' sim_global = plspm(simdata, sim_inner, 
#'                    sim_outer, modes=sim_mod)
#' sim_global
#' 
#' ## Then compute cluster analysis on residuals of global model
#' sim_clus = res.clus(sim_global)
#' 
#' ## To complete REBUS, run iterative algorithm
#' rebus_sim = it.reb(sim_global, sim_clus, nk=2, 
#'                    stop.crit=0.005, iter.max=100)
#' 
#' ## You can also compute complete outputs 
#' ## for local models by running:
#' local_rebus = local.models(sim_global, rebus_sim)
#' 
#' # Display plspm summary for first local model
#' summary(local_rebus$loc.model.1)
#' }
it.reb <- 
function(pls, hclus.res, nk, Y = NULL, stop.crit = 0.005, iter.max = 100)
{
  # =======================================================
  # checking arguments
  # =======================================================
  if (class(pls) != "plspm") 
    stop("\n'it.reb()' requires an object of class 'plspm'")
  if (any(pls$model$specs$modes != "A"))
    stop("\nREBUS only works for reflective modes")
  if (!pls$model$specs$scaled)
    stop("\nREBUS only works with scaled='TRUE'")
  if (missing(hclus.res))
    stop("\nargument 'hclus.res' is missing")
  if (class(hclus.res) != "hclust")
    stop("\n'it.reb()' requires an object of class 'hclust'")
  if (missing(nk))
    stop("\nargument 'nk' (number of classes) is missing")
  if (mode(nk) != "numeric" || length(nk) != 1 || 
        nk <= 1 || (nk%%1) != 0)
    stop("\nInvalid number of classes 'nk'. Must be an integer larger than 1")
  # test availibility of dataset (either Y or pls$data)
  test_dataset(Y, pls$data, pls$model$gens$obs)
  # stop.crit
  if (mode(stop.crit) != "numeric" || length(stop.crit) != 1 || 
        stop.crit < 0 || stop.crit >= 1)
  {
    warning("Invalid stop criterion 'stop.crit'. Deafult value 0.005 is used")
    stop.crit = 0.005
  }
  if (mode(iter.max) != "numeric" || length(iter.max) != 1 || 
        iter.max <= 1 || (iter.max%%1) != 0)
  {
    warning("Invalid 'iter.max'. Deafult value 100 is used")
    iter.max = 100
  }
  
  # =======================================================
  # inputs settings
  # =======================================================
  IDM = pls$model$IDM
  blocks = pls$model$blocks
  specs = pls$model$specs
  blocklist = indexify(blocks)
  # data matrix DM
  if (!is.null(pls$data)) {
    DM = pls$data
    if (!is.matrix(DM)) DM = as.matrix(DM)
    dataset = TRUE
  } else {         
    dataset = FALSE
    # building data matrix 'DM'
    DM = get_manifests(Y, blocks)
  }
  lvs = pls$model$gens$lvs
  lvs.names = pls$model$gens$lvs_names
  mvs = pls$model$gens$mvs
  mvs.names = pls$model$gens$mvs_names
  N = pls$model$gens$obs
  endo = rowSums(IDM)
  endo[endo != 0] = 1  
  n.end = sum(endo) 
  # data scaling
  X = get_data_scaled(DM, TRUE)
  # initial partition (cutting dendrogram in 'nk' clusters)
  ini.part = cutree(hclus.res, nk)
  # number of clusters
  nclus = nlevels(factor(ini.part)) 
  
  # =======================================================
  # initialize variables to store results
  # =======================================================
  w.locals <- as.list(1:nclus)# outer.weights
  Y.locals <- as.list(1:nclus)# std latent variables 
  loads.locals <- as.list(1:nclus)# loadings
  path.locals <- as.list(1:nclus)# path coefficients
  R2.locals <- as.list(1:nclus)# R2
  comu.locals <- as.list(1:nclus)# mvs communalities
  outres.locals <- as.list(1:nclus)# communality residuals
  innres.locals <- as.list(1:nclus)# structural residuals
  CM.locals <- matrix(NA, N, nclus)# CM distance of each unit from each model
  
  # =======================================================
  # iterative process
  # =======================================================
  old.clas = ini.part
  iter.ch = 0
  repeat 
  {
    # define MV matrix for each initial class
    split.DM <- split.X <- as.list(1:nclus)
    for (k in 1:nclus) split.DM[[k]] = DM[old.clas==k,]            
    # local models computation
    for (k in 1:nclus)
    {   
      nk = nrow(split.DM[[k]])
      # local mean
      mean.k = apply(split.DM[[k]], 2, mean)
      # local std.dev
      sd.k = sqrt((nk-1)/nk) * apply(split.DM[[k]], 2, sd)
      # spliting data matrix for each class
      split.X[[k]] = scale(split.DM[[k]], center=mean.k, scale=sd.k)
      # calculating outer weights for each class
      out.ws = get_weights(split.X[[k]], IDM, blocks, specs)
      w.locals[[k]] = out.ws$W
      # calculating LV scores for each class
      Y.k = split.X[[k]] %*% out.ws$W
      # calculating path coefficients for each class
      pathmod = get_paths(IDM, Y.k)
      path.locals[[k]] = pathmod[[2]]
      R2.locals[[k]] = pathmod[[3]][endo==1]
      # calculating loadings and communalities for each class
      xloads = cor(split.X[[k]], Y.k)
      loads.locals[[k]] = rowSums(xloads * out.ws$ODM)
      comu.locals[[k]] = loads.locals[[k]]^2
      
      # latent variables for each unit in each local model
      X.k = scale(DM, center=mean.k, scale=sd.k)
      Y.locals[[k]] = X.k %*% w.locals[[k]]
      # computation of communality residuals
      out.res = DM
      for (j in 1:lvs)
      {
        X.hat = Y.locals[[k]][,j] %*% t(loads.locals[[k]][blocklist==j])
        # outer residuals
        out.res[,blocklist==j] = (X.k[,blocklist==j] - X.hat)^2
      }
      outres.locals[[k]] = out.res
      # computation of inner residuals
      if (sum(endo) != 1)
        Y.hat = Y.locals[[k]] %*% t(path.locals[[k]][endo==1,])   
      if (sum(endo) == 1)
        Y.hat = Y.locals[[k]] %*% path.locals[[k]][endo==1,]        
      innres.locals[[k]] = (Y.locals[[k]][,endo==1] - Y.hat)^2
      # computation super normalized residual of outer model
      res.num1 = outres.locals[[k]] %*% diag(1/comu.locals[[k]], mvs, mvs)
      supres.outer = rowSums(res.num1) / (sum(rowSums(res.num1))/(N-2))
      # computation of super normalized residual of inner model
      res.num2 = innres.locals[[k]] %*% diag(1/R2.locals[[k]],n.end,n.end)
      supres.inner = rowSums(res.num2) / (sum(rowSums(res.num2))/(N-2)) 
      # computation of the CM distance
      CM.locals[,k] = sqrt(supres.outer * supres.inner)
    }
    # allocating the units to their closest class
    new.clas = old.clas
    for (i in 1:N)
      new.clas[i] <- which(CM.locals[i,] == min(CM.locals[i,]))[1]
    # checking convergence
    dif.clas = new.clas - old.clas 
    unit.change = length(which(dif.clas != 0))
    iter.ch = iter.ch + 1
    # rate of unit change
    if (unit.change/N < stop.crit || iter.ch == iter.max)
      break
    old.clas = new.clas
    if (any(table(new.clas) <= 5)) 
      stop("\nToo few units: a class with less than 6 units was detected") 
  }
  
  # =======================================================
  # computation of final local models
  # =======================================================
  # data matrices for each class
  DM.cl.rebus = as.list(1:nclus)    
  for (k in 1:nclus)
    DM.cl.rebus[[k]] = DM[new.clas==k,]
  # data matrix with units membership
  DM.rebus = cbind(DM, rebus.class=new.clas)
  redun.locals = as.list(1:nclus)
  # table of path coefficients for local models
  path.labs = NULL
  for (j in 1:lvs)
    for (i in j:lvs)
      if (IDM[i,j] == 1) 
        path.labs = c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
  reb.labs = paste(rep("Class",nclus), 1:nclus, sep=".")
  path.rebus = matrix(NA, length(path.labs), nclus)
  loads.rebus = matrix(NA, mvs, nclus)
  qual.rebus = matrix(NA, sum(lvs+2*n.end+1), nclus)
  for (k in 1:nclus)
  {               
    nk = nrow(DM.cl.rebus[[k]])
    # local mean
    mean.k = apply(DM.cl.rebus[[k]], 2, mean)
    # local std.dev
    sd.k = sqrt((nk-1)/nk) * apply(DM.cl.rebus[[k]], 2, sd)
    DM.k = scale(DM.cl.rebus[[k]], center=mean.k, scale=sd.k)        
    out.ws = get_weights(DM.k, IDM, blocks, specs)
    Y.k = DM.k %*% out.ws$W
    pathmod = get_paths(IDM, Y.k)
    path.locals[[k]] = pathmod[[2]]
    R2.locals[[k]] = pathmod[[3]][endo==1]
    xloads = cor(DM.cl.rebus[[k]], Y.k)
    loads.locals[[k]] = rowSums(xloads * out.ws$ODM)
    comu.locals[[k]] = loads.locals[[k]]^2
    redun = rep(0, mvs)
    aux = 0
    for (j in 1:lvs) {
      if (endo[j] == 1) {
        aux = aux + 1
        redun_aux = comu.locals[[k]][blocklist==j] * R2.locals[[k]][aux]
        redun[blocklist==j] = redun_aux
      }
    }
    redun.locals[[k]] = redun   
    # path coeffs for local models
    path.rebus[,k] = as.vector(path.locals[[k]][IDM==1])
    # loadings for local models
    loads.rebus[,k] = loads.locals[[k]]
    # table of quality indexes
    comu.aveg = rep(NA, lvs) 
    redun.aveg <- R2.aux <- rep(NA, n.end)
    aux = 0
    for (j in 1:lvs)
    {
      comu.aveg[j] = mean(comu.locals[[k]][blocklist==j]) 
      if (endo[j] == 1) 
      {
        aux = aux + 1      
        redun.aveg[aux] = mean(redun.locals[[k]][blocklist==j])
        R2.aux[aux] = R2.locals[[k]][aux]
      }
    }
    qual.rebus[1:lvs,k] = comu.aveg
    qual.rebus[(lvs+1):(lvs+n.end),k] = redun.aveg
    qual.rebus[(lvs+n.end+1):(lvs+2*n.end),k] = R2.aux
    qual.rebus[(lvs+2*n.end+1),k] = sqrt(mean(comu.aveg) * mean(R2.aux))
  }
  gqi = get_GQI(pls, new.clas, DM)
  dimnames(path.rebus) = list(path.labs, reb.labs)
  dimnames(loads.rebus) = list(mvs.names, reb.labs)
  v1 = paste(rep("Com", lvs), lvs.names, sep=".")
  v2 = paste(rep("Red", n.end), lvs.names[endo==1], sep=".")
  v3 = paste(rep("R2", n.end), lvs.names[endo==1], sep=".")
  dimnames(qual.rebus) = list(c(v1,v2,v3,"GoF"), reb.labs)
  qual.rebus = qual.rebus
  aux = list(lvs, n.end, unit.change/N, stop.crit, iter.max, iter.ch, gqi)
  res = list(path.coef = path.rebus, 
             loadings = loads.rebus, 
             quality = qual.rebus,
             segments = new.clas, 
             origdata.clas = DM.rebus, 
             aux = aux)
  class(res) = "rebus"
  return(res)
}
