\name{get_paths}
\alias{get_paths}
\title{Calculate path coefficients for \code{plspm}}
\usage{
  get_paths(IDM, Y.lvs, plsr)
}
\arguments{
  \item{IDM}{Inner Design Matrix}

  \item{Y.lvs}{Matrix of latent variables}

  \item{plsr}{logical to indicate path coefs by pls
  regression}
}
\value{
  list with inner mode, path coefs matrix, R2, and
  residuals
}
\description{
  Internal function. \code{get_paths} is called by
  \code{plspm}.
}
\keyword{internal}

