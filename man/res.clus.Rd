\name{res.clus}
\alias{res.clus}
\title{Cluster Analysis on communality and structural
residuals of a PLS-PM object}
\usage{
  res.clus(pls, Y = NULL)
}
\arguments{
  \item{pls}{Object of class \code{"plspm"}}

  \item{Y}{Optional dataset (matrix or data frame) used
  when argument \code{dataset=NULL} inside \code{pls}.}
}
\value{
  An Object of class \code{"hclust"} containing the results
  of the Hierarchical Cluster Analysis on the communality
  and structural residuals.
}
\description{
  Computes communality and structural residuals from the
  global model and performs a Hierarchical Cluster Analysis
  on these residuals according to the REBUS algorithm.
}
\details{
  The function \code{res.clus} comprises the second and
  third steps of the REBUS-PLS Algorithm. It computes
  communality and structural residuals. Then it performs a
  Hierarchical Cluster Analysis on these residuals (step
  three of REBUS-PLS Algorithm). As a result, this function
  directly provides a dendrogram obtained from a
  Hierarchical Cluster Analysis.
}
\examples{
\dontrun{
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

 # Then compute cluster analysis on the residuals of global model
 sim_clus = res.clus(sim_global)
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
}
\seealso{
  \code{\link{it.reb}}, \code{\link{plspm}}
}

