\name{local.models}
\alias{local.models}
\title{PLS-PM for global and local models}
\usage{
  local.models(pls, y, Y = NULL)
}
\arguments{
  \item{pls}{An object of class \code{"plspm"}}

  \item{y}{One object of the following classes:
  \code{"rebus"}, \code{"integer"}, or \code{"factor"},
  that provides the class partitions.}

  \item{Y}{Optional dataset (matrix or data frame) used
  when argument \code{dataset=NULL} inside \code{pls}.}
}
\value{
  An object of class \code{"local.models"}, basically a
  list of length \code{k+1}, where \code{k} is the number
  of classes.

  \item{glob.model}{PLS-PM of the global model}

  \item{loc.model.1}{PLS-PM of segment (class) 1}

  \item{loc.model.2}{PLS-PM of segment (class) 2}

  \item{loc.model.k}{PLS-PM of segment (class) k}
}
\description{
  Calculates PLS-PM for global and local models from a
  given partition
}
\details{
  \code{local.models} calculates PLS-PM for the global
  model (i.e. over all observations) as well as PLS-PM for
  local models (i.e. observations of different partitions).

  When \code{y} is an object of class \code{"rebus"},
  \code{local.models} is applied to the classes obtained
  from the REBUS algorithm.

  When \code{y} is an \code{integer} vector or a
  \code{factor}, the values or levels are assumed to
  represent the group to which each observation belongs. In
  this case, the function \code{local.models} calculates
  PLS-PM for the global model, as well as PLS-PM for each
  group (local models).

  When the object \code{pls} does not contain a data matrix
  (i.e. \code{pls$data=NULL}), the user must provide the
  data matrix or data frame in \code{Y}.

  The original parameters \code{modes}, \code{scheme},
  \code{scaled}, \code{tol}, and \code{iter} from the
  object \code{pls} are taken.
}
\note{
  Each element of the list is an object of class
  \code{"plspm"}. Thus, in order to examine the results for
  each local model, it is necessary to use the
  \code{summary} function.
}
\examples{
## Example of REBUS PLS with simulated data

# load simdata
data(simdata, package="plspm2")

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
\author{
  Laura Trinchera, Gaston Sanchez
}
\seealso{
  \code{\link{rebus.pls}}
}

