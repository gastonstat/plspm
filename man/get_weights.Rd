\name{get_weights}
\alias{get_weights}
\title{Calculate outer weights for \code{plspm}}
\usage{
  get_weights(X, IDM, blocks, modes, scheme, tol, iter)
}
\arguments{
  \item{X}{scaled data}

  \item{IDM}{Inner Design Matrix}

  \item{blocks}{list with variables in each block}

  \item{modes}{vector measurement mode}

  \item{scheme}{inner path scheme to use}

  \item{tol}{tolerance threshold for convergen}

  \item{iter}{maximum number of iterations}
}
\description{
  Internal function. \code{get_weights} is called by
  \code{plspm}
}
\keyword{internal}

