\name{get_boots}
\alias{get_boots}
\title{Performs bootstrap validation in \code{plspm}}
\usage{
  get_boots(DM, path_matrix, blocks, specs, br)
}
\arguments{
  \item{DM}{Data Matrix}

  \item{path_matrix}{Inner Design Matrix}

  \item{blocks}{list of vectors with column indices}

  \item{specs}{list with algorithm specifications}

  \item{br}{number of bootstrap resamples}
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{get_boots} is called by
  \code{plspm}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

