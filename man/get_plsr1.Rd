\name{get_plsr1}
\alias{get_plsr1}
\title{PLS regression for \code{plspm}}
\usage{
  get_plsr1(predictors, response, nc = NULL, scaled = TRUE)
}
\arguments{
  \item{predictors}{matrix of predictors}

  \item{response}{response variable}

  \item{nc}{number of components}

  \item{scaled}{logical indicating whether to scale the
  data}
}
\value{
  A list with pls regression results
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_plsr1} performs PLS-R1.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

