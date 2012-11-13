\name{plspm}
\alias{plspm}
\title{PLS-PM: Partial Least Squares Path Modeling}
\usage{
  plspm(Data, inner_matrix, outer_list, modes = NULL,
    scheme = "centroid", scaled = TRUE, tol = 1e-05,
    iter = 100, boot.val = FALSE, br = NULL, plsr = FALSE,
    dataset = TRUE)
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

  \item{boot.val}{A logical value indicating whether
  bootstrap validation is performed (\code{FALSE} by
  default).}

  \item{br}{An integer indicating the number bootstrap
  resamples. Used only when \code{boot.val=TRUE}. When
  \code{boot.val=TRUE}, the default number of re-samples is
  100, but it can be specified in a range from 100 to
  1000.}

  \item{plsr}{A logical value indicating whether pls
  regression is applied to calculate path coefficients
  (\code{FALSE} by default).}

  \item{dataset}{A logical value indicating whether the
  data matrix should be retrieved (\code{TRUE} by
  default).}
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

  \item{outer.cor}{Correlations between the latent
  variables and the manifest variables (also called
  crossloadings)}

  \item{inner.sum}{Summarized results by latent variable of
  the inner model. Includes: type of LV, type of
  measurement, number of indicators, R-squared, average
  communality, average redundancy, and average variance
  extracted}

  \item{effects}{Path effects of the structural
  relationships. Includes: direct, indirect, and total
  effects}

  \item{unidim}{Results for checking the unidimensionality
  of blocks (These results are only meaningful for
  reflective blocks)}

  \item{gof}{Goodness-of-Fit index}

  \item{data}{Data matrix containing the manifest variables
  used in the model. Only when \code{dataset=TRUE}}

  \item{boot}{List of bootstrapping results; only available
  when argument \code{boot.val=TRUE}}
}
\description{
  Estimate path models with latent variables by partial
  least squares approach
}
\details{
  The function \code{plspm} estimates a path model by
  partial least squares approach providing the full set of
  results. \cr

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

 # vector of modes (reflective indicators)
 sat_mod = rep("A", 6)

 # apply plspm
 satpls = plspm(satisfaction, sat_inner, sat_outer, sat_mod, scaled=FALSE, boot.val=TRUE)

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

  Tenenhaus M., Pages J. (2001) Multiple factor analysis
  combined with PLS path modelling. Application to the
  analysis of relationships between physicochemical
  variables, sensory profiles and hedonic judgements.
  \emph{Chemometrics and Intelligent Laboratory Systems},
  \bold{58}, pp. 261-273.

  Tenenhaus M., Hanafi M. (2010) A bridge between PLS path
  modeling and multi-block data analysis. \emph{Handbook on
  Partial Least Squares (PLS): Concepts, methods, and
  applications.} Springer.

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
\seealso{
  \code{\link{plspm.fit}}, \code{\link{plot.plspm}}
}

