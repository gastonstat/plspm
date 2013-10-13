\name{get_nom_scale}
\alias{get_nom_scale}
\title{Non-Metric Nominal Scale}
\usage{
  get_nom_scale(y, x, Xdummy)
}
\arguments{
  \item{y}{vector of values}

  \item{x}{vector of the natural number series 1, 2, ..., p
  obtained by means of the function myRank}

  \item{Xdummy}{dummy matrix corresponding to x}
}
\value{
  scaled matrix
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_nom_scale} is called by
  \code{plspm}.
}
\note{
  This function replaces the elements of x by the the means
  of y conditioned to the levels of x
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

