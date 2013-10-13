\name{get_gof}
\alias{get_gof}
\title{Goodness-of-fit for \code{plspm}}
\usage{
  get_gof(comu, R2, blocks, path_matrix)
}
\arguments{
  \item{comu}{list of communalities}

  \item{R2}{vector of R-squared coefficients}

  \item{blocks}{list of variables in each block}

  \item{path_matrix}{Inner Design Matrix}
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_gof} is called by
  \code{plspm}
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

