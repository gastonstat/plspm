\name{get_PLSR_NA}
\alias{get_PLSR_NA}
\title{Internal PLS regression with missing values}
\usage{
  get_PLSR_NA(Y, X, ncomp)
}
\arguments{
  \item{Y}{an (already centered) vector of order n;
  \code{NA} not allowed}

  \item{X}{an (already centered) n x p matrix; \code{NA}
  allowed but with at least one observation for every row
  and column}

  \item{ncomp}{number of components}
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_PLSR_NA}
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

