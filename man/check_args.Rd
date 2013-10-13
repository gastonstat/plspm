\name{check_args}
\alias{check_args}
\title{Check arguments for \code{plspm} and \code{plspm.fit}}
\usage{
  check_args(Data, path_matrix, blocks, scaling, modes,
    scheme, scaled, tol, maxiter, plscomp, boot.val, br,
    dataset)
}
\arguments{
  \item{Data}{numeric matrix or data frame containing the
  manifest variables.}

  \item{path_matrix}{square (lower triangular) boolean
  matrix for inner model.}

  \item{blocks}{List of vectors (either numeric or
  character) to specify the blocks of manifest variables
  asociated to the latent variables.}

  \item{modes}{character indicating the type of
  measurement.}

  \item{scheme}{string indicating the type of inner
  weighting scheme.}

  \item{scaled}{logical indicating whether scaling data is
  performed.}

  \item{tol}{decimal value indicating the tol criterion for
  covergence.}

  \item{maxiter}{integer indicating the maximum number of
  iterations.}

  \item{plscomp}{optional vector indicating the number of
  PLS components (for each block) to be used when handling
  non-metric data (only used if \code{scaling} is
  provided)}

  \item{bootstrap}{whether bootstrap validation is
  performed.}

  \item{br}{integer indicating the number bootstrap
  resamples.}

  \item{dataset}{whether the data matrix should be
  retrieved.}
}
\value{
  list of validated arguments
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_args} is called by
  \code{plspm}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

