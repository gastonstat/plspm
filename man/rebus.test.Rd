\name{rebus.test}
\alias{rebus.test}
\title{Permutation Test for REBUS Multi-Group Comparison}
\usage{
  rebus.test(pls, reb, Y = NULL)
}
\arguments{
  \item{pls}{Object of class \code{"plspm"} returned by
  \code{\link{plspm}}}

  \item{reb}{Object of class \code{"rebus"} returned by
  either \code{\link{rebus.pls}} or \code{\link{it.reb}}.}

  \item{Y}{Optional dataset (matrix or data frame) used
  when argument \code{dataset=NULL} inside \code{pls}.}
}
\value{
  An object of class \code{"rebus.test"}, basically a list
  containing the results of each pair of compared classes.
  In turn, each element of the list is also a list with the
  results for the path coefficients, loadings, and GoF
  index.
}
\description{
  Performs permutation tests for comparing pairs of groups
  from a REBUS object.
}
\details{
  A permutation test on path coefficients, loadings, and
  GoF index is applied to the classes obtained from REBUS,
  by comparing two classes at a time. That is to say, a
  permutation test is applied on pair of classes. The
  number of permutations in each test is 100. In turn, the
  number of classes handled by \code{rebus.test} is limited
  to 6.

  When \code{pls$data=NULL} (there is no data matrix), the
  user must provide the data matrix or data frame in
  \code{Y}.
}
\examples{
\dontrun{
 ## typical example of PLS-PM in customer satisfaction analysis
 ## model with six LVs and reflective indicators
 ## example of rebus analysis with simulated data

 # load data
 data(simdata)

 # Calculate plspm
 sim_path = matrix(c(0,0,0,0,0,0,1,1,0), 3, 3, byrow=TRUE)
 dimnames(sim_path) = list(c("Price", "Quality", "Satisfaction"),
                            c("Price", "Quality", "Satisfaction"))
 sim_blocks = list(c(1,2,3,4,5), c(6,7,8,9,10), c(11,12,13))
 sim_mod = c("A", "A", "A")  # reflective indicators
 sim_global = plspm(simdata, sim_path,
                    sim_blocks, modes=sim_mod)
 sim_global

 # Cluster analysis on residuals of global model
 sim_clus = res.clus(sim_global)

 # Iterative steps of REBUS algorithm on 2 classes
 rebus_sim = it.reb(sim_global, sim_clus, nk=2,
                    stop.crit=0.005, iter.max=100)

 # apply rebus.test
 sim_permu = rebus.test(sim_global, rebus_sim)

 # inspect sim.rebus
 sim_permu
 sim_permu$test_1_2

 # or equivalently
 sim_permu[[1]]
 }
}
\author{
  Laura Trinchera, Gaston Sanchez
}
\references{
  Chin, W.W. (2003) A permutation procedure for multi-group
  comparison of PLS models. In: Vilares M., Tenenhaus M.,
  Coelho P., Esposito Vinzi V., Morineau A. (Eds.)
  \emph{PLS and Related Methods - Proceedings of the
  International Symposium PLS03.} Decisia, pp. 33-43.
}
\seealso{
  \code{\link{rebus.pls}}, \code{\link{local.models}}
}

