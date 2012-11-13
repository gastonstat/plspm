\name{get_unidim}
\alias{get_unidim}
\title{Unidimensionality of reflective blocks}
\usage{
  get_unidim(x = NULL, outer = NULL, modes = NULL,
    DM = NULL, blocks = NULL, check = TRUE)
}
\arguments{
  \item{x}{numeric matrix or data frame containing the
  manifest variables.}

  \item{outer}{list of vectors with column indices from
  \code{x} indicating the sets of manifest variables
  asociated to the latent variables}

  \item{modes}{character vector indicating the type of
  measurement}

  \item{DM}{Data Matrix}

  \item{blocks}{vector with numbers of variables per block}

  \item{check}{logical to check \code{x}, \code{outer} and
  \code{modes}.}
}
\value{
  A data frame with the following columns:

  \item{Type.measure}{Measurement mode}

  \item{MVs}{number of manifest variables in each block}

  \item{C.alpha}{Cronbach's alpha}

  \item{DG.rho}{Dillon-Goldstein rho}

  \item{eig.1st}{First eigenvalue}

  \item{eig.2nd}{Second eigenvalue}
}
\description{
  Internal this function. \code{get_unidim} is called by
  \code{plspm}
}
\keyword{internal}

