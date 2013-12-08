\name{it.reb}
\alias{it.reb}
\title{Iterative steps of Response-Based Unit Segmentation (REBUS)}
\usage{
  it.reb(pls, hclus.res, nk, Y = NULL, stop.crit = 0.005,
    iter.max = 100)
}
\arguments{
  \item{pls}{an object of class \code{"plspm"}}

  \item{hclus.res}{object of class \code{"res.clus"}
  returned by \code{\link{res.clus}}}

  \item{nk}{integer larger than 1 indicating the number of
  classes. This value should be defined according to the
  dendrogram obtained by performing
  \code{\link{res.clus}}.}

  \item{Y}{optional data matrix used when \code{pls$data}
  is \code{NULL}}

  \item{stop.crit}{Number indicating the stop criterion for
  the iterative algorithm. It is suggested to use the
  threshold of less than 0.05\% of units changing class
  from one iteration to the other as stopping rule.}

  \item{iter.max}{integer indicating the maximum number of
  iterations}
}
\value{
  an object of class \code{"rebus"}

  \item{loadings}{Matrix of standardized loadings (i.e.
  correlations with LVs.) for each local model}

  \item{path.coefs}{Matrix of path coefficients for each
  local model}

  \item{quality}{Matrix containing the average
  communalities, the average redundancies, the R2 values,
  and the GoF index for each local model}

  \item{segments}{Vector defining the class membership of
  each unit}

  \item{origdata.clas}{ The numeric matrix with original
  data and with a new column defining class membership of
  each unit}
}
\description{
  REBUS-PLS is an iterative algorithm for performing
  response based clustering in a PLS-PM framework.
  \code{it.reb} allows to perform the iterative steps of
  the REBUS-PLS Algorithm. It provides summarized results
  for final local models and the final partition of the
  units. Before running this function, it is necessary to
  run the \code{\link{res.clus}} function to choose the
  number of classes to take into account.
}
\examples{
\dontrun{
## Example of REBUS PLS with simulated data

# load simdata
data("simdata", package='plspm')

# Calculate global plspm
sim_inner = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
dimnames(sim_inner) = list(c("Price", "Quality", "Satisfaction"),
                           c("Price", "Quality", "Satisfaction"))
sim_outer = list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13))
sim_mod = c("A", "A", "A")  # reflective indicators
sim_global = plspm(simdata, sim_inner,
                   sim_outer, modes=sim_mod)
sim_global

## Then compute cluster analysis on residuals of global model
sim_clus = res.clus(sim_global)

## To complete REBUS, run iterative algorithm
rebus_sim = it.reb(sim_global, sim_clus, nk=2,
                   stop.crit=0.005, iter.max=100)

## You can also compute complete outputs
## for local models by running:
local_rebus = local.models(sim_global, rebus_sim)

# Display plspm summary for first local model
summary(local_rebus$loc.model.1)
}
}
\author{
  Laura Trinchera, Gaston Sanchez
}
\references{
  Esposito Vinzi, V., Trinchera, L., Squillacciotti, S.,
  and Tenenhaus, M. (2008) REBUS-PLS: A Response-Based
  Procedure for detecting Unit Segments in PLS Path
  Modeling. \emph{Applied Stochastic Models in Business and
  Industry (ASMBI)}, \bold{24}, pp. 439-458.

  Trinchera, L. (2007) Unobserved Heterogeneity in
  Structural Equation Models: a new approach to latent
  class detection in PLS Path Modeling. \emph{Ph.D.
  Thesis}, University of Naples "Federico II", Naples,
  Italy.

  \url{http://www.fedoa.unina.it/view/people/Trinchera,_Laura.html}
}
\seealso{
  \code{\link{plspm}}, \code{\link{rebus.pls}},
  \code{\link{res.clus}}
}

