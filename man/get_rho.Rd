\name{get_rho}
\alias{get_rho}
\title{Calculate Dillon-Goldstein's rho}
\usage{
  get_rho(Block, score)
}
\arguments{
  \item{Block}{matrix (one block) of standardized manifest
  variables}

  \item{score}{first principal component of \code{Block}}
}
\value{
  Dillon-Goldstein's rho
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_rho} is called by
  \code{get_unidim} and \code{rho}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

