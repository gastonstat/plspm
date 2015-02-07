\docType{package}
\name{plspm}
\alias{plspm}
\alias{plspm-package}
\title{PLS-PM: Partial Least Squares Path Modeling}
\usage{
  plspm(Data, path_matrix, blocks, modes = NULL,
    scaling = NULL, scheme = "centroid", scaled = TRUE,
    tol = 1e-06, maxiter = 100, plscomp = NULL,
    boot.val = FALSE, br = NULL, dataset = TRUE)
}
\arguments{
  \item{Data}{matrix or data frame containing the manifest
  variables.}

  \item{path_matrix}{A square (lower triangular) boolean
  matrix representing the inner model (i.e. the path
  relationships between latent variables).}

  \item{blocks}{list of vectors with column indices or
  column names from \code{Data} indicating the sets of
  manifest variables forming each block (i.e. which
  manifest variables correspond to each block).}

  \item{scaling}{optional list of string vectors indicating
  the type of measurement scale for each manifest variable
  specified in \code{blocks}. \code{scaling} must be
  specified when working with non-metric variables.
  Possible values: \code{"num"} (numeric), \code{"raw"},
  \code{"nom"} (nominal), and \code{"ord"} (ordinal).}

  \item{modes}{character vector indicating the type of
  measurement for each block. Possible values are:
  \code{"A", "B", "newA", "PLScore", "PLScow"}. The length
  of \code{modes} must be equal to the length of
  \code{blocks}.}

  \item{scheme}{string indicating the type of inner
  weighting scheme. Possible values are \code{"centroid"},
  \code{"factorial"}, or \code{"path"}.}

  \item{scaled}{whether manifest variables should be
  standardized. Only used when \code{scaling = NULL}. When
  (\code{TRUE}, data is scaled to standardized values
  (mean=0 and variance=1). The variance is calculated
  dividing by \code{N} instead of \code{N-1}).}

  \item{tol}{decimal value indicating the tolerance
  criterion for the iterations (\code{tol=0.000001}). Can
  be specified between 0 and 0.001.}

  \item{maxiter}{integer indicating the maximum number of
  iterations (\code{maxiter=100} by default). The minimum
  value of \code{maxiter} is 100.}

  \item{plscomp}{optional vector indicating the number of
  PLS components (for each block) to be used when handling
  non-metric data (only used if \code{scaling} is
  provided)}

  \item{boot.val}{whether bootstrap validation should be
  performed. (\code{FALSE} by default).}

  \item{br}{number bootstrap resamples. Used only when
  \code{boot.val=TRUE}. When \code{boot.val=TRUE}, the
  default number of re-samples is 100, but it can be
  specified in a range from 100 to 1000.}

  \item{dataset}{whether the data matrix used in the
  computations should be retrieved (\code{TRUE} by
  default).}
}
\value{
  An object of class \code{"plspm"}.

  \item{outer_model}{Results of the outer model. Includes:
  outer weights, standardized loadings, communalities, and
  redundancies}

  \item{inner_model}{Results of the inner (structural)
  model. Includes: path coeffs and R-squared for each
  endogenous latent variable}

  \item{scores}{Matrix of latent variables used to estimate
  the inner model. If \code{scaled=FALSE} then
  \code{scores} are latent variables calculated with the
  original data (non-stardardized).}

  \item{path_coefs}{Matrix of path coefficients (this
  matrix has a similar form as \code{path_matrix})}

  \item{crossloadings}{Correlations between the latent
  variables and the manifest variables (also called
  crossloadings)}

  \item{inner_summary}{Summarized results of the inner
  model. Includes: type of LV, type of measurement, number
  of indicators, R-squared, average communality, average
  redundancy, and average variance extracted}

  \item{effects}{Path effects of the structural
  relationships. Includes: direct, indirect, and total
  effects}

  \item{unidim}{Results for checking the unidimensionality
  of blocks (These results are only meaningful for
  reflective blocks)}

  \item{gof}{Goodness-of-Fit index}

  \item{data}{Data matrix containing the manifest variables
  used in the model. Only available when
  \code{dataset=TRUE}}

  \item{boot}{List of bootstrapping results; only available
  when argument \code{boot.val=TRUE}}
}
\description{
  Estimate path models with latent variables by partial
  least squares approach (for both metric and non-metric
  data)

  Estimate path models with latent variables by partial
  least squares approach (for both metric and non-metric
  data)
}
\details{
  The function \code{plspm} estimates a path model by
  partial least squares approach providing the full set of
  results. \cr

  The argument \code{path_matrix} is a matrix of zeros and
  ones that indicates the structural relationships between
  latent variables. \code{path_matrix} must be a lower
  triangular matrix; it contains a 1 when column \code{j}
  affects row \code{i}, 0 otherwise. \cr

  \itemize{ \item \code{\link{plspm}}: Partial Least
  Squares Path Modeling \item \code{\link{plspm.fit}}:
  Simple version for PLS-PM \item
  \code{\link{plspm.groups}}: Two Groups Comparison in
  PLS-PM \item \code{\link{rebus.pls}}: Response Based Unit
  Segmentation (REBUS) }
}
\examples{
\dontrun{
## typical example of PLS-PM in customer satisfaction analysis
## model with six LVs and reflective indicators

# load dataset satisfaction
data(satisfaction)

# path matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
sat_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

# plot diagram of path matrix
innerplot(sat_path)

# blocks of outer model
sat_blocks = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

# vector of modes (reflective indicators)
sat_mod = rep("A", 6)

# apply plspm
satpls = plspm(satisfaction, sat_path, sat_blocks, modes = sat_mod,
   scaled = FALSE)

# plot diagram of the inner model
innerplot(satpls)

# plot loadings
outerplot(satpls, what = "loadings")

# plot outer weights
outerplot(satpls, what = "weights")
}
}
\author{
  Gaston Sanchez, Giorgio Russolillo
}
\references{
  Tenenhaus M., Esposito Vinzi V., Chatelin Y.M., and Lauro
  C. (2005) PLS path modeling. \emph{Computational
  Statistics & Data Analysis}, \bold{48}, pp. 159-205.

  Lohmoller J.-B. (1989) \emph{Latent variables path
  modeling with partial least squares.} Heidelberg:
  Physica-Verlag.

  Wold H. (1985) Partial Least Squares. In: Kotz, S.,
  Johnson, N.L. (Eds.), \emph{Encyclopedia of Statistical
  Sciences}, Vol. 6. Wiley, New York, pp. 581-591.

  Wold H. (1982) Soft modeling: the basic design and some
  extensions. In: K.G. Joreskog & H. Wold (Eds.),
  \emph{Systems under indirect observations: Causality,
  structure, prediction}, Part 2, pp. 1-54. Amsterdam:
  Holland.

  Russolillo, G. (2012) Non-Metric Partial Least Squares.
  \emph{Electronic Journal of Statistics}, \bold{6}, pp.
  1641-1669.
  \url{http://projecteuclid.org/euclid.ejs/1348665231}
}
\seealso{
  \code{\link{innerplot}}, \code{\link{outerplot}},
}

