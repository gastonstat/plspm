\name{check_modes}
\alias{check_modes}
\title{Check modes}
\usage{
  check_modes(modes, blocks)
}
\arguments{
  \item{modes}{character vector indicating the type of
  measurement for each block. Values regular user can
  choose: \code{"A", "B", "newA"} Values Expert user can
  choose: \code{"A", "B", "newA", "PLScore", "PLScow"}. The
  length of \code{modes} must be equal to the length of
  \code{blocks}.}

  \item{blocks}{list defining the blocks of manifest
  variables}
}
\value{
  validated modes
}
\description{
  Internal function not to be called by the user
}
\details{
  Internal function. \code{check_modes} is called by
  \code{check_specs}.
}
\section{Warning}{
  Do NOT use this function unless you are ME, a package
  developer, or a jedi user who really knows what is doing
  (seriously!)
}
\keyword{internal}

