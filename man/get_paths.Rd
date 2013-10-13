\name{get_paths}
\alias{get_paths}
\title{Calculate path coefficients for \code{plspm}}
\usage{
  get_paths(path_matrix, Y_lvs, full = TRUE)
}
\arguments{
  \item{path_matrix}{path matrix}

  \item{Y_lvs}{Matrix of latent variables}

  \item{full}{logical to indicate all results from
  'summary(lm())'}
}
\value{
  list with inner results, path coefs matrix, R2, and
  residuals
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_paths} is called by
  \code{plspm}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

