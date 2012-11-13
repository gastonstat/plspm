\name{rebus.pls}
\alias{rebus.pls}
\title{Response Based Unit Segmentation (REBUS)}
\usage{
  rebus.pls(pls, Y = NULL, stop.crit = 0.005,
    iter.max = 100)
}
\arguments{
  \item{pls}{Object of class \code{"plspm"}}

  \item{Y}{Optional dataset (matrix or data frame) used
  when argument \code{dataset=NULL} inside \code{pls}.}

  \item{stop.crit}{Number indicating the stop criterion for
  the iterative algorithm. Use a threshold of less than
  0.05\% of units changing class from one iteration to the
  other as stopping rule.}

  \item{iter.max}{integer indicating the maximum number of
  iterations.}
}
\value{
  An object of class \code{"rebus"}, basically a list with:

  \item{loadings}{Matrix of standardized loadings (i.e.
  correlations with LVs.) for each local model.}

  \item{path.coefs}{Matrix of path coefficients for each
  local model.}

  \item{quality}{Matrix containing the average
  communalities, average redundancies, R2 values, and GoF
  values for each local model.}

  \item{segments}{Vector defining for each unit the class
  membership.}

  \item{origdata.clas}{The numeric matrix with original
  data and with a new column defining class membership of
  each unit.}
}
\description{
  Performs all the steps of the REBUS-PLS algorithm.
  Starting from the global model, REBUS allows us to detect
  local models with better performance.
}
\examples{
\dontrun{
 ## typical example of PLS-PM in customer satisfaction analysis
 ## model with six LVs and reflective indicators
 ## example of rebus analysis with simulated data

 # load data
 data(simdata)

 # Calculate plspm
 sim_inner = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
 dimnames(sim_inner) = list(c("Price", "Quality", "Satisfaction"),
                            c("Price", "Quality", "Satisfaction"))
 sim_outer = list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13))
 sim_mod = c("A", "A", "A")  # reflective indicators
 sim_global = plspm(simdata, inner=sim_inner,
                    outer=sim_outer, modes=sim_mod)
 sim_global

 # run rebus.pls and choose the number of classes
 # to be taken into account according to the displayed dendrogram.
 rebus_sim = rebus.pls(sim_global, stop.crit = 0.005, iter.max = 100)

 # You can also compute complete outputs for local models by running:
 local_rebus = local.models(sim_global, rebus_sim)
 }
}
\author{
  Laura Trinchera, Gaston Sanchez
}
\references{
  Esposito Vinzi V., Trinchera L., Squillacciotti S., and
  Tenenhaus M. (2008) REBUS-PLS: A Response-Based Procedure
  for detecting Unit Segments in PLS Path Modeling.
  \emph{Applied Stochastic Models in Business and Industry
  (ASMBI)}, \bold{24}, pp. 439-458.

  Trinchera, L. (2007) Unobserved Heterogeneity in
  Structural Equation Models: a new approach to latent
  class detection in PLS Path Modeling. \emph{Ph.D.
  Thesis}, University of Naples "Federico II", Naples,
  Italy.

  \url{http://www.fedoa.unina.it/view/people/Trinchera,_Laura.html}
}
\seealso{
  \code{\link{plspm}}, \code{\link{res.clus}},
  \code{\link{it.reb}}, \code{\link{rebus.test}},
  \code{\link{local.models}}
}

