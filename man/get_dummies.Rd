\name{get_dummies}
\alias{get_dummies}
\title{Dummy matrices for categorical manifest variables}
\usage{
  get_dummies(MV, specs)
}
\arguments{
  \item{MV}{matrix or data frame from with manifest
  variables}

  \item{specs}{list with algorithm specifications}
}
\value{
  dummies list with dummy matrices for categorical manifest
  variables
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_dummies} is called by
  \code{plspm}.
}
\note{
  If there are no categorical variables specified in specs,
  \code{get_dummies} returns \code{NULL}
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

