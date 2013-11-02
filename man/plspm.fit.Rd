\name{plspm.fit}
\alias{plspm.fit}
\title{Basic results for Partial Least Squares Path Modeling}
\usage{
  plspm.fit(Data, path_matrix, blocks, modes = NULL,
    scaling = NULL, scheme = "centroid", scaled = TRUE,
    tol = 1e-06, maxiter = 100, plscomp = NULL)
}
\arguments{
  \item{Data}{matrix or data frame containing the manifest
  variables.}

  \item{path_matrix}{A square (lower triangular) boolean
  matrix representing the inner model (i.e. the path
  relationships betwenn latent variables).}

  \item{blocks}{list of vectors with column indices or
  column names from \code{Data} indicating the sets of
  manifest variables forming each block (i.e. which
  manifest variables correspond to each block).}

  \item{scaling}{optional list of string vectors indicating
  the type of measurement scale for each manifest variable
  specified in \code{blocks}. \code{scaling} must be
  specified when working with non-metric variables.}

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
  (\code{TRUE} data is scaled to standardized values
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
  original data (non-stardardized). If \code{scaled=TRUE}
  then \code{scores} and \code{latents} have the same
  values}

  \item{path_coefs}{Matrix of path coefficients (this
  matrix has a similar form as \code{path_matrix})}
}
\description{
  Estimate path models with latent variables by partial
  least squares approach without providing the full list of
  results as \code{plspm()}. This might be helpful when
  doing simulations, intensive computations, or when you
  don't want the whole enchilada.
}
\details{
  \code{plspm.fit} performs the basic PLS algorithm and
  provides limited results (e.g. outer model, inner model,
  scores, and path coefficients). \cr

  The argument \code{path_matrix} is a matrix of zeros and
  ones that indicates the structural relationships between
  latent variables. \code{path_matrix} must be a lower
  triangular matrix; it contains a 1 when column \code{j}
  affects row \code{i}, 0 otherwise. \cr
}
\examples{
\dontrun{
 ## typical example of PLS-PM in customer satisfaction analysis
 ## model with six LVs and reflective indicators

 # load dataset satisfaction
 data(satisfaction)

 # inner model matrix
 IMAG = c(0,0,0,0,0,0)
 EXPE = c(1,0,0,0,0,0)
 QUAL = c(0,1,0,0,0,0)
 VAL = c(0,1,1,0,0,0)
 SAT = c(1,1,1,1,0,0)
 LOY = c(1,0,0,0,1,0)
 sat_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

 # outer model list
 sat_blocks = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

 # vector of reflective modes
 sat_modes = rep("A", 6)

 # apply plspm.fit
 satpls = plspm.fit(satisfaction, sat_path, sat_blocks, sat_modes,
     scaled=FALSE)

 # summary of results
 summary(satpls)

 # default plot (inner model)
 plot(satpls)
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
}
\seealso{
  \code{\link{innerplot}}, \code{\link{plot.plspm}},
}

