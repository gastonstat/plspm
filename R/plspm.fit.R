#'@title Basic results for Partial Least Squares Path Modeling
#'
#'@description
#'Estimate path models with latent variables by partial least squares approach
#'without providing the full list of results as \code{plspm}. This might be helpful
#'when doing simulations, intensive computations, or when you don't want
#'the whole enchilada.
#'
#'@details
#'\code{plspm.fit} performs the basic PLS algorithm and provides
#'limited results (e.g. outer weights, LVs scores, path coefficients, R2, and
#'loadings). \cr
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
#'@author Gaston Sanchez
#'
#'@references Tenenhaus M., Esposito Vinzi V., Chatelin Y.M., and Lauro C.
#'(2005) PLS path modeling. \emph{Computational Statistics & Data Analysis},
#'\bold{48}, pp. 159-205.
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
#'  # vector of reflective modes
#'  sat_mod = rep("A", 6)
#'
#'  # apply plspm.fit
#'  satpls = plspm.fit(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE)
#'  
#'  # summary of results
#'  summary(satpls)
#'
#'  # default plot (inner model)
#'  plot(satpls)
#'  }
#'
plspm.fit <-
function(Data, inner_matrix, outer_list, modes = NULL, scheme="centroid", 
         scaled = TRUE, tol = 0.00001, iter = 100)
{
    # =======================================================
    # checking arguments
    # =======================================================
    valid = get_params(x=Data, inner=inner_matrix, outer=outer_list, modes=modes, 
                       scheme=scheme, scaled=scaled, boot.val=FALSE, br=NULL, 
                       plsr=FALSE, tol=tol, iter=iter, dataset=FALSE)
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
    # inner design matrix
    IDM = inner
    # latent variables names
    lvs.names = rownames(IDM)
    dimnames(IDM) = list(lvs.names, lvs.names)
    # how many latent variables
    lvs = nrow(IDM)
    # get blocks: number of mvs per block
    blocks = unlist(lapply(outer, length))
    mvs = sum(blocks)
    names(blocks) = lvs.names
    blocklist = outer
    for (k in 1:length(outer))
         blocklist[[k]] <- rep(k,blocks[k])
    blocklist = unlist(blocklist)
    Mode = modes
    Mode[modes=="A"] = "Reflective"
    Mode[modes=="B"] = "Formative"   
    # building data matrix 'DM'
    DM = matrix(NA, nrow(x), mvs)
    mvs.names = rep(NA, mvs)
    for (k in 1:lvs)
    {        
        DM[,which(blocklist==k)] = as.matrix(x[,outer[[k]]])
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
    out.ws <- get_weights(X, IDM, blocks, modes, scheme, tol, iter)
    if (is.null(out.ws)) {
        print(paste("Iterative process is non-convergent with 'iter'=", 
                    iter, " and 'tol'=", tol, sep=""))
        stop("Algorithm stops") 
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
    # =======================================================
    # Stage 3: Measurement loadings and communalities
    # =======================================================  
    loadcomu = get_loads(X, Y.lvs, blocks)    
    loads = loadcomu[[1]]
    comu = loadcomu[[2]]

    # =======================================================
    # Measurement model
    # =======================================================  
    outmod <- as.list(1:lvs)
    for (j in 1:lvs)
    {
        aux <- which(blocklist==j)
        outmod[[j]] <- round(cbind(weights=out.weights[aux], std.loads=loads[aux], 
                             communal=comu[aux]), 4)
    }
    names(outmod) = lvs.names

    # =======================================================
    # Results
    # =======================================================  
    skem = switch(scheme, "centroid"="centroid", "factor"="factor", "path"="path")
    model = list(IDM=IDM, blocks=blocks, scheme=skem, modes=modes, scaled=scaled, 
                  obs=nrow(X), tol=tol, iter=iter, n.iter=out.ws[[3]], outer=outer)
    res = list(outer.mod = outmod, 
               inner.mod = innmod, 
               latents = Z.lvs, 
               scores = Y.lvs,
               out.weights = out.weights, 
               loadings = loads, 
               path.coefs = Path, 
               r.sqr = R2, 
               data = NULL, 
               model = model)
    class(res) = c("plspm.fit", "plspm")
    return(res)
}

