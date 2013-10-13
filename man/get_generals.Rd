\name{get_generals}
\alias{get_generals}
\title{Get general parameters}
\usage{
  get_generals(MV, path_matrix)
}
\arguments{
  \item{MV}{matrix of manifest variables}

  \item{path_matrix}{matrix with path connections}
}
\value{
  list with number and names of observations, MVs and LVs
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_generals} is called by
  \code{plspm}. This function gets the number and names of:
  observations, manifest variables, and latent variables
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

