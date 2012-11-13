\name{get_params}
\alias{get_params}
\title{Check parameters for \code{plspm} and \code{plspm.fit}}
\usage{
  get_params(x, inner, outer, modes = NULL,
    scheme = "centroid", scaled = TRUE, boot.val = FALSE,
    br = NULL, plsr = FALSE, tol = 1e-05, iter = 100,
    dataset = TRUE)
}
\arguments{
  \item{x}{numeric matrix or data frame containing the
  manifest variables.}

  \item{inner}{square (lower triangular) boolean matrix for
  inner model.}

  \item{outer}{List of vectors with column indices from
  \code{x} indicating the sets of manifest variables
  asociated to the latent variables.}

  \item{modes}{character vector indicating the type of
  measurement.}

  \item{scheme}{string indicating the type of inner
  weighting scheme.}

  \item{scaled}{logical indicating whether scaling data is
  performed.}

  \item{boot.val}{logical indicating whether bootstrap
  validation is performed.}

  \item{br}{integer indicating the number bootstrap
  resamples.}

  \item{plsr}{logical indicating whether to use pls
  regression for path coefs.}

  \item{tol}{decimal value indicating the tolerance
  criterion for covergence.}

  \item{iter}{integer indicating the maximum number of
  iterations.}

  \item{dataset}{logical indicating whether the data matrix
  should be retrieved.}
}
\value{
  list of validated parameters for \code{plspm} and
  \code{plspm.fit}.
}
\description{
  Internal function. \code{get_params} is called by
  \code{plspm}.
}
\keyword{internal}

