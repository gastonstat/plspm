#'@title PLS-PM: Partial Least Squares Path Modeling
#'
#'@description
#'Estimate path models with latent variables by partial least squares approach
#'
#'@details
#'The function \code{plspm} estimates a path model by partial least squares
#'approach providing the full set of results. \cr
#'
#'The argument \code{inner_matrix} is a matrix of zeros and ones that indicates
#'the structural relationships between latent variables. This must be a lower
#'triangular matrix. \code{inner_matrix} will contain a 1 when column \code{j}
#'affects row \code{i}, 0 otherwise. \cr
#'
#'@param Data A numeric matrix or data frame containing the manifest variables.
#'@param inner_matrix A square (lower triangular) boolean matrix representing the
#'inner model (i.e. the path relationships betwenn latent variables).
#'@param outer_list List of vectors with column indices from \code{x} indicating 
#'the sets of manifest variables asociated to the latent variables
#'(i.e. which manifest variables correspond to the latent variables).
#'Length of \code{outer_list} must be equal to the number of rows of \code{inner_matrix}.
#'@param modes A character vector indicating the type of measurement for each
#'latent variable. \code{"A"} for reflective measurement or \code{"B"} for
#'formative measurement (\code{NULL} by default). The length of \code{modes}
#'must be equal to the length of \code{outer_list}).
#'@param scheme A string of characters indicating the type of inner weighting
#'scheme. Possible values are \code{"centroid"}, \code{"factor"}, or
#'\code{"path"}.
#'@param scaled A logical value indicating whether scaling data is performed
#'When (\code{TRUE} data is scaled to standardized values (mean=0 and variance=1)
#'The variance is calculated dividing by \code{N} instead of \code{N-1}).
#'@param tol Decimal value indicating the tolerance criterion for the
#'iterations (\code{tol=0.00001}). Can be specified between 0 and 0.001.
#'@param iter An integer indicating the maximum number of iterations
#'(\code{iter=100} by default). The minimum value of \code{iter} is 100.
#'@param boot.val A logical value indicating whether bootstrap validation is
#'performed (\code{FALSE} by default). 
#'@param br An integer indicating the number bootstrap resamples. Used only
#'when \code{boot.val=TRUE}. When \code{boot.val=TRUE}, the default number of 
#'re-samples is 100, but it can be specified in a range from 100 to 1000.
#'@param plsr A logical value indicating whether pls regression is applied
#'to calculate path coefficients (\code{FALSE} by default).
#'@param dataset A logical value indicating whether the data matrix should be
#'retrieved (\code{TRUE} by default).
#'@return An object of class \code{"plspm"}. 
#'@return \item{outer.mod}{Results of the outer (measurement) model. Includes:
#'outer weights, standardized loadings, communalities, and redundancies}
#'@return \item{inner.mod}{Results of the inner (structural) model. Includes: path
#'coefficients and R-squared for each endogenous latent variable}
#'@return \item{latents}{Matrix of standardized latent variables (variance=1
#'calculated divided by \code{N}) obtained from centered data (mean=0)}
#'@return \item{scores}{Matrix of latent variables used to estimate the inner
#'model. If \code{scaled=FALSE} then \code{scores} are latent variables
#'calculated with the original data (non-stardardized). If \code{scaled=TRUE}
#'then \code{scores} and \code{latents} have the same values}
#'@return \item{out.weights}{Vector of outer weights}
#'@return \item{loadings}{Vector of standardized loadings (i.e. correlations with
#'LVs)}
#'@return \item{path.coefs}{Matrix of path coefficients (this matrix has a similar
#'form as \code{inner_matrix})}
#'@return \item{r.sqr}{Vector of R-squared coefficients}
#'@return \item{outer.cor}{Correlations between the latent variables and the
#'manifest variables (also called crossloadings)}
#'@return \item{inner.sum}{Summarized results by latent variable of the inner
#'model. Includes: type of LV, type of measurement, number of indicators,
#'R-squared, average communality, average redundancy, and average variance
#'extracted}
#'@return \item{effects}{Path effects of the structural relationships. Includes:
#'direct, indirect, and total effects}
#'@return \item{unidim}{Results for checking the unidimensionality of blocks
#'(These results are only meaningful for reflective blocks)}
#'@return \item{gof}{Goodness-of-Fit index}
#'@return \item{data}{Data matrix containing the manifest variables used in the
#'model. Only when \code{dataset=TRUE}}
#'@return \item{boot}{List of bootstrapping results; only available when argument
#'\code{boot.val=TRUE}}
#'@author Gaston Sanchez
#'
#'@references Tenenhaus M., Esposito Vinzi V., Chatelin Y.M., and Lauro C.
#'(2005) PLS path modeling. \emph{Computational Statistics & Data Analysis},
#'\bold{48}, pp. 159-205.
#'
#'Tenenhaus M., Pages J. (2001) Multiple factor analysis combined with
#'PLS path modelling. Application to the analysis of relationships between
#'physicochemical variables, sensory profiles and hedonic judgements.
#'\emph{Chemometrics and Intelligent Laboratory Systems}, \bold{58}, pp.
#'261-273.
#'
#'Tenenhaus M., Hanafi M. (2010) A bridge between PLS path modeling and
#'multi-block data analysis. \emph{Handbook on Partial Least Squares (PLS):
#'Concepts, methods, and applications.} Springer.
#'
#'Lohmoller J.-B. (1989) \emph{Latent variables path modelin with partial
#'least squares.} Heidelberg: Physica-Verlag.
#'
#'Wold H. (1985) Partial Least Squares. In: Kotz, S., Johnson, N.L. (Eds.),
#'\emph{Encyclopedia of Statistical Sciences}, Vol. 6. Wiley, New York, pp.
#'581-591.
#'
#'Wold H. (1982) Soft modeling: the basic design and some extensions. In: K.G.
#'Joreskog & H. Wold (Eds.), \emph{Systems under indirect observations:
#'Causality, structure, prediction}, Part 2, pp. 1-54. Amsterdam: Holland.
#'@seealso \code{\link{plspm.fit}}, \code{\link{plot.plspm}}
#'@export
#'@examples
#'
#'  \dontrun{
#'  ## typical example of PLS-PM in customer satisfaction analysis
#'  ## model with six LVs and reflective indicators
#'
#'  # load dataset satisfaction
#'  data(satisfaction)
#'
#'  # inner model matrix
#'  IMAG = c(0,0,0,0,0,0)
#'  EXPE = c(1,0,0,0,0,0)
#'  QUAL = c(0,1,0,0,0,0)
#'  VAL = c(0,1,1,0,0,0)
#'  SAT = c(1,1,1,1,0,0) 
#'  LOY = c(1,0,0,0,1,0)
#'  sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
#'
#'  # outer model list
#'  sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
#'
#'  # vector of modes (reflective indicators)
#'  sat_mod = rep("A", 6)
#'
#'  # apply plspm
#'  satpls = plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=TRUE)
#'  
#'  # summary of results
#'  summary(satpls)
#'
#'  # default plot (inner model)
#'  plot(satpls)
#'  }
#'
plspm <-
function(Data, inner_matrix, outer_list, modes = NULL, scheme = "centroid", 
         scaled = TRUE, tol = 0.00001, iter = 100, boot.val = FALSE, 
         br = NULL, plsr = FALSE, dataset = TRUE)
{
  # =======================================================
  # checking arguments
  # =======================================================
  valid = get_params(x=Data, inner=inner_matrix, outer=outer_list, modes=modes, 
                     scheme=scheme, scaled=scaled, tol=tol, iter=iter,
                     boot.val=boot.val, br=br, plsr=plsr, dataset=dataset)
  x = valid$x
  inner = valid$inner
  outer = valid$outer
  modes = valid$modes
  scheme = valid$scheme
  SCHEMES = valid$SCHEMES
  scaled = valid$scaled
  boot.val = valid$boot.val
  br = valid$br
  plsr = valid$plsr
  tol = valid$tol
  iter = valid$iter
  dataset = valid$dataset
  
  # =======================================================
  # inputs setting
  # =======================================================  
  IDM = inner
  lvs.names = rownames(IDM)
  dimnames(IDM) = list(lvs.names, lvs.names)
  lvs = nrow(IDM)
  blocks = unlist(lapply(outer, length))
  mvs = sum(blocks)
  names(blocks) = lvs.names
  blocklist = outer
  for (k in 1:length(blocks))
    blocklist[[k]] = rep(k,blocks[k])
  blocklist = unlist(blocklist)
  Mode = modes
  Mode[modes=="A"] = "Reflective"
  Mode[modes=="B"] = "Formative"   
  # building data matrix 'DM'
  DM = matrix(NA, nrow(x), mvs)
  mvs.names = rep(NA, mvs)
  for (k in 1:lvs)
  {        
    DM[,which(blocklist==k)] <- as.matrix(x[,outer[[k]]])
    mvs.names[which(blocklist==k)] = colnames(x)[outer[[k]]]
  }
  dimnames(DM) = list(rownames(x), mvs.names)
  # apply the selected scaling
  if (scaled) {
    sd.X = sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
    X = scale(DM, scale=sd.X)
  } else {
    X = scale(DM, scale=FALSE)
  }
  dimnames(X) = list(rownames(x), mvs.names)
  
  # =======================================================
  # Stage 1: Iterative procedure
  # =======================================================  
  out.ws = get_weights(X, IDM, blocks, modes, scheme, tol, iter)
  if (is.null(out.ws)) {
    print(paste("Iterative process is non-convergent with 'iter'=", 
                iter, " and 'tol'=", tol, sep=""))
    message("Algorithm stops") 
    stop("")
  }
  out.weights = out.ws[[1]]
  cor.XY = cor(X, X%*%out.ws[[2]])
  w.sig = rep(NA, lvs)
  for (k in 1:lvs) 
    w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
  Z.lvs = X %*% out.ws[[2]] %*% diag(w.sig, lvs, lvs)
  Y.lvs = Z.lvs
  if (!scaled) 
    Y.lvs = DM %*% out.ws[[2]] %*% diag(w.sig, lvs, lvs)
  dimnames(Y.lvs) = list(rownames(X), lvs.names)
  dimnames(Z.lvs) = list(rownames(X), lvs.names)
  # =======================================================
  # Stage 2: Path coefficients and total effects
  # =======================================================  
  pathmod = get_paths(IDM, Y.lvs, plsr)
  innmod = pathmod[[1]]
  Path = pathmod[[2]]
  R2 = pathmod[[3]]
  Path.efs = get_effects(Path)
  # =======================================================
  # Stage 3: Measurement loadings and communalities
  # =======================================================  
  loadcomu = get_loads(X, Y.lvs, blocks)    
  loads = loadcomu[[1]]
  comu = loadcomu[[2]]
  endo = rowSums(IDM)
  endo[endo != 0] = 1  
  redun = rep(0, mvs)
  for (j in 1:lvs)
    if (endo[j] == 1)
      redun[blocklist==j] = comu[blocklist==j] * R2[j]
  
  # =======================================================
  # Measurement model
  # =======================================================  
  outcor <- outmod <- as.list(1:lvs)
  for (j in 1:lvs)
  {
    aux <- which(blocklist==j)
    outmod[[j]] <- round(cbind(weights=out.weights[aux], std.loads=loads[aux], 
                               communal=comu[aux], redundan=redun[aux]), 4)
    outcor[[j]] <- round(cor(DM[,aux], Y.lvs), 4)
  }
  names(outmod) = lvs.names
  names(outcor) = lvs.names  
  
  # =======================================================
  # Unidimensionality
  # ======================================================= 
  unidim = get_unidim(x=NULL, outer=NULL, modes=modes, 
                         DM=DM, blocks=blocks, check=FALSE)
  
  # =======================================================
  # Summary Inner model
  # =======================================================  
  exo.endo = rowSums(IDM)
  exo.endo[rowSums(IDM)==0] = "Exogen"
  exo.endo[rowSums(IDM)!=0] = "Endogen"
  av.comu = rep(0, lvs)   # average communality
  av.redu = rep(0, lvs)   # average redundancy
  ave = rep(0, lvs)      # average variance extracted
  for (k in 1:lvs)
  {
    av.comu[k] = mean(comu[which(blocklist==k)])
    av.redu[k] = mean(redun[which(blocklist==k)])
    if (modes[k]=="A")
    {
      ave.num = sum(comu[which(blocklist==k)])
      ave.denom = sum(comu[which(blocklist==k)]) + sum(1-(comu[which(blocklist==k)]))
      ave[k] = ave.num / ave.denom
    }
  }
  names(ave) = lvs.names
  innsum = data.frame(LV.Type = exo.endo, 
                      Measure = abbreviate(Mode, 5), 
                      MVs = blocks, 
                      R.square = R2, 
                      Av.Commu = av.comu, 
                      Av.Redun = av.redu, 
                      AVE = ave)
  rownames(innsum) = lvs.names
  
  # =======================================================
  # GoF Index
  # =======================================================  
  gof = get_gof(comu, R2, blocks, IDM)
  
  # =======================================================
  # Results
  # =======================================================  
  skem <- switch(scheme, "centroid"="centroid", "factor"="factor", "path"="path")
  model <- list(IDM=IDM, blocks=blocks, scheme=skem, modes=modes, scaled=scaled, 
                boot.val=boot.val, plsr=plsr, obs=nrow(X), br=br, 
                tol=tol, iter=iter, n.iter=out.ws[[3]], outer=outer)
  # deliver dataset?
  if (dataset) data = DM else data = NULL
  # deliver bootstrap validation results? 
  if (boot.val) 
  {
    if (nrow(X) <= 10) {
      warning("Bootstrapping stopped: very few cases.") 
    } else { 
      n.efs = nrow(Path.efs)
      res.boot = get_boots(DM, IDM, blocks, modes, scheme, scaled, br, plsr, tol, iter)
    }
  } else {
    res.boot = FALSE
  }
  res = list(outer.mod = outmod, 
              inner.mod = innmod, 
              latents = Z.lvs, 
              scores = Y.lvs,
              out.weights = out.weights, 
              loadings = loads, 
              path.coefs = Path, 
              r.sqr = R2,
              outer.cor = outcor, 
              inner.sum = innsum, 
              effects = Path.efs,
              unidim = unidim, 
              gof = gof, 
              boot = res.boot, 
              data = data, 
              model = model)
  class(res) = "plspm"
  return(res)
}





