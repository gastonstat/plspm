\name{check_tol}
\alias{check_tol}
\title{Check tolerance threshold}
\usage{
  check_tol(tol)
}
\arguments{
  \item{tol}{decimal value indicating the tolerance
  threshold for the convergence of the PLS algorithm
  (iterative stage)}
}
\value{
  validated tol
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_tol} is called by
  \code{check_specs}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

