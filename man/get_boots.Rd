\name{get_boots}
\alias{get_boots}
\title{Performs bootstrap validation in \code{plspm}}
\usage{
  get_boots(DM, IDM, blocks, modes, scheme, scaled, br,
    plsr, tol, iter)
}
\arguments{
  \item{DM}{Data Matrix}

  \item{IDM}{Inner Design Matrix}

  \item{blocks}{vector with number of lvs per block}

  \item{modes}{vector of modes}

  \item{scheme}{type of inner path scheme}

  \item{scaled}{logical to scale the data}

  \item{br}{number of bootstrap resamples}

  \item{plsr}{logical to get path coefficients by pls
  regression}

  \item{tol}{threshold for convergence}

  \item{iter}{maximum number of iterations}
}
\description{
  Internal function. \code{get_boots} is called by
  \code{plspm}.
}
\keyword{internal}

