\name{check_specs}
\alias{check_specs}
\title{Check specifications of PLS-PM algorithm}
\usage{
  check_specs(blocks, scaling, modes, scheme, scaled, tol,
    maxiter, plscomp)
}
\arguments{
  \item{blocks}{list defining the blocks of manifest
  variables}

  \item{scaling}{list with type of scaling for each
  manifest variable}

  \item{modes}{character vector specifying the modes}

  \item{scheme}{character string indicating the inner
  weighting scheme}

  \item{scaled}{should manifest variables be standardized}

  \item{tol}{character indicating the inner weighting
  scheme}

  \item{maxiter}{character indicating the inner weighting
  scheme}

  \item{plscomp}{optional vector indicating the number of
  PLS components (for each block) to be used when handling
  non-metric data (only used if \code{scaling} is
  provided)}
}
\value{
  validated specifications
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_specs} is called by
  \code{check_args}.
}
\note{
  if scaling is specified, it overrides scaled
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

