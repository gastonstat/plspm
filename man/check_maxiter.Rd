\name{check_maxiter}
\alias{check_maxiter}
\title{Check maximum number of iterations}
\usage{
  check_maxiter(maxiter)
}
\arguments{
  \item{maxiter}{integer indicating the maximum number of
  iterations (iterative stage of the PLS algorithm)}
}
\value{
  validated maxiter
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_maxiter} is called by
  \code{check_specs}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

