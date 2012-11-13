\name{plspm.fit}
\alias{plspm.fit}
\title{Basic results for Partial Least Squares Path Modeling}
\usage{
  plspm.fit(Data, inner_matrix, outer_list, modes = NULL,
    scheme = "centroid", scaled = TRUE, tol = 1e-05,
    iter = 100)
}
\arguments{
  \item{Data}{A numeric matrix or data frame containing the
  manifest variables.}

  \item{inner_matrix}{A square (lower triangular) boolean
  matrix representing the inner model (i.e. the path
  relationships betwenn latent variables).}

  \item{outer_list}{List of vectors with column indices
  from \code{x} indicating the sets of manifest variables
  asociated to the latent variables (i.e. which manifest
  variables correspond to the latent variables). Length of
  \code{outer_list} must be equal to the number of rows of
  \code{inner_matrix}.}

  \item{modes}{A character vector indicating the type of
  measurement for each latent variable. \code{"A"} for
  reflective measurement or \code{"B"} for formative
  measurement (\code{NULL} by default). The length of
  \code{modes} must be equal to the length of
  \code{outer_list}).}

  \item{scheme}{A string of characters indicating the type
  of inner weighting scheme. Possible values are
  \code{"centroid"}, \code{"factor"}, or \code{"path"}.}

  \item{scaled}{A logical value indicating whether scaling
  data is performed When (\code{TRUE} data is scaled to
  standardized values (mean=0 and variance=1) The variance
  is calculated dividing by \code{N} instead of
  \code{N-1}).}

  \item{tol}{Decimal value indicating the tolerance
  criterion for the iterations (\code{tol=0.00001}). Can be
  specified between 0 and 0.001.}

  \item{iter}{An integer indicating the maximum number of
  iterations (\code{iter=100} by default). The minimum
  value of \code{iter} is 100.}
}
\value{
  An object of class \code{"plspm"}.

  \item{outer.mod}{Results of the outer (measurement)
  model. Includes: outer weights, standardized loadings,
  communalities, and redundancies}

  \item{inner.mod}{Results of the inner (structural) model.
  Includes: path coefficients and R-squared for each
  endogenous latent variable}

  \item{latents}{Matrix of standardized latent variables
  (variance=1 calculated divided by \code{N}) obtained from
  centered data (mean=0)}

  \item{scores}{Matrix of latent variables used to estimate
  the inner model. If \code{scaled=FALSE} then
  \code{scores} are latent variables calculated with the
  original data (non-stardardized). If \code{scaled=TRUE}
  then \code{scores} and \code{latents} have the same
  values}

  \item{out.weights}{Vector of outer weights}

  \item{loadings}{Vector of standardized loadings (i.e.
  correlations with LVs)}

  \item{path.coefs}{Matrix of path coefficients (this
  matrix has a similar form as \code{inner_matrix})}

  \item{r.sqr}{Vector of R-squared coefficients}
}
\description{
  Estimate path models with latent variables by partial
  least squares approach without providing the full list of
  results as \code{plspm}. This might be helpful when doing
  simulations, intensive computations, or when you don't
  want the whole enchilada.
}
\details{
  \code{plspm.fit} performs the basic PLS algorithm and
  provides limited results (e.g. outer weights, LVs scores,
  path coefficients, R2, and loadings). \cr

  The argument \code{inner_matrix} is a matrix of zeros and
  ones that indicates the structural relationships between
  latent variables. This must be a lower triangular matrix.
  \code{inner_matrix} will contain a 1 when column \code{j}
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
 sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)

 # outer model list
 sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

 # vector of reflective modes
 sat_mod = rep("A", 6)

 # apply plspm.fit
 satpls = plspm.fit(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE)

 # summary of results
 summary(satpls)

 # default plot (inner model)
 plot(satpls)
 }
}
\author{
  Gaston Sanchez
}
\references{
  Tenenhaus M., Esposito Vinzi V., Chatelin Y.M., and Lauro
  C. (2005) PLS path modeling. \emph{Computational
  Statistics & Data Analysis}, \bold{48}, pp. 159-205.

  Lohmoller J.-B. (1989) \emph{Latent variables path
  modelin with partial least squares.} Heidelberg:
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

